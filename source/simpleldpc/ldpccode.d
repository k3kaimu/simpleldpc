module simpleldpc.ldpccode;

import std.algorithm : min, max;
import std.exception : enforce;
import std.format : format;
import std.math : tanh, atanh;
import std.range : iota;

import simpleldpc.wifildpc;

/++
BLDPC encoder and decoder.
+/
class BLDPCCode
{
    this(uint block_length, uint info_length)
    {
        _N = block_length;
        _K = info_length;
        _Z = 0;

        _M = _N - _K;
        _H_mat = new ubyte[][](_N, _M);
    }


    /**
    パリティ検査行列Hとcodewordをかけ合わせて，パリティ検査します
    */
    bool checkCodeword(in ubyte[] decoded_cw) const
    {
        bool check = true;
        foreach(i_check; 0 .. _M) {
            ubyte c;
            foreach(i_col; _row_mat[i_check])
                c += decoded_cw[i_col];
            
            if(c % 2 == 1) {
                check = false;
                break;
            }
        }

        return check;
    }


    void loadWiFiLDPC(uint block_length, uint rate_index)
    {
        enum BLOCK_LENS = [648, 1296, 1944];
        enum ZS = [27, 54, 81];
        LswitchBlockLen1: switch(block_length)
        {
          static foreach(i, BLOCK_LEN; BLOCK_LENS)
          {
            case BLOCK_LEN:
                _Z = ZS[i];
                break LswitchBlockLen1;
          }

            default:
                enforce(false, "Block length value not supported for WiFi LDPC");
                return;
        }

        _N = block_length;

        enum RATES = [[1, 2], [2, 3], [3, 4], [5, 6]];

        LswitchRateIndex: switch(rate_index)
        {
            static foreach(i, RATE; RATES)
            {
                case i:
                    _K = _N * RATE[0] / RATE[1];

                    mixin(q{LswitchBlockLen2_%1$s: switch(block_length)
                    {
                        static foreach(j, BLOCK_LEN; BLOCK_LENS)
                        {
                            case BLOCK_LEN:
                                generateQCLDPC(mixin("WiFiLDPC.H_%%s_%%s_%%s".format(BLOCK_LEN, RATE[0], RATE[1])));
                                break LswitchBlockLen2_%1$s;
                        }

                        default:
                            enforce(false, "Block length value not supported for WiFi LDPC");
                            return;
                    }}.format(i));
                    break LswitchRateIndex;
            }

            default:
                enforce(false, "Rate index value not supported");
                return;
        }

        _M = _N - _K;
    }


    uint infoLength() const
    {
        return _K;
    }


    /++
    S. Myung, K. Yang, J. Kim, "Quasi-Cyclic LDPC Codes for Fast Encoding", IEEE Trans. on Inform. Theory, vol. 51, no. 8, August 2005.
    +/
    ubyte[] encode(in ubyte[] info_bits) const
    {
        ubyte[] codeword = new ubyte[](_N);
        codeword[0 .. _K] = info_bits[];

        ubyte[] parity = new ubyte[_M];

        // H_I u^T
        foreach(i_row, row; _row_mat) {
            foreach(i_col; row) if(i_col < _K)
                parity[i_row] += codeword[i_col];
            
            parity[i_row] %= 2;
        }

        // p_1 = ET^{-1} As + Cs
        foreach(i_col; 0 .. _Z) {
            foreach(i_row; iota(i_col, _M, _Z))
                codeword[_K + i_col] += parity[i_row];

            codeword[_K + i_col] %= 2;
        }

        // compute As + Bp_1
        foreach(i_row, row; _row_mat) {
            foreach(i_col; row) if(i_col >= _K && i_col < _K + _Z)
                parity[i_row] += codeword[i_col];

            parity[i_row] %= 2;
        }

        // back substitution
        foreach(i_col; iota(_K + _Z, _N, _Z)) {
            foreach(i_row; 0 .. _Z) {
                codeword[i_col + i_row] = parity[i_col + i_row - _K - _Z];
                parity[i_col + i_row - _K] += parity[i_col + i_row - _K - _Z];
            }
        }

        foreach(ref e; codeword) e %= 2;

        return codeword;
    }


    ubyte[] decodeLLR(in double[] llr, uint max_iter) const
    {
        double[][] edge_mat = new double[][](_M);
        double[][] last_edge_mat = new double[][](_M);
        double[] updated_llr = llr.dup;

        ubyte[] decoded_cw = new ubyte[_N];

        foreach(i; 0 .. _M) {
            edge_mat[i] = new double[](_row_mat[i].length);
            edge_mat[i][] = 0;
            last_edge_mat[i] = new double[](_row_mat[i].length);
            last_edge_mat[i][] = 0;
        }

        foreach(iter; 0 .. max_iter) 
        {
            // 関数ノードでの周辺化
            foreach(i_row, row; _row_mat) {
                // i_rowのチェックノードからi_col_1の変数ノードへ送るLLRを計算する
                foreach(i_col_index1, i_col_1; row) {
                    double tmp = 1;
                    foreach(i_col_index2, i_col_2; row) {
                        if(i_col_index1 == i_col_index2) continue;

                        // 変数ノードでの計算時に自身のノードのLLRも加えているので
                        // それは引き算して周辺化した値を得る
                        double l1 = updated_llr[i_col_2] - last_edge_mat[i_row][i_col_index2];
                        // tanhの発散を防ぐ
                        l1 = min(l1, 20.0);
                        l1 = max(l1, -20.0);
                        tmp *= tanh(l1/2);
                    }

                    edge_mat[i_row][i_col_index1] = 2 * atanh(tmp);
                }
            }

            // copy
            foreach(i; 0 .. _M)
                last_edge_mat[i][] = edge_mat[i][];

            // 変数ノードでの計算では，事前分布としての入力LLRに関数ノードからのLLRを加える
            // そのためここではまず入力LLRを入れておく
            updated_llr[] = llr[];

            // 関数ノードからのメッセージの合成積，ただし，対数尤度比なので和を計算する
            // 自分のノード分も足しているので，それは関数ノードでの計算で引き算する
            foreach(i_row, row; _row_mat) {
                foreach(i_col_index, i_col; row) {
                    updated_llr[i_col] += last_edge_mat[i_row][i_col_index];
                }
            }

            // LLRに従って硬判定をする
            foreach(i; 0 .. _N) {
                if(updated_llr[i] > 0)
                    decoded_cw[i] = 0;
                else
                    decoded_cw[i] = 1;
            }

            if(checkCodeword(decoded_cw))
                break;
        }

        return decoded_cw;
    }

  private:
    ubyte[][] _H_mat;       // transposed
    uint _N;
    uint _K;
    uint _M;
    uint _Z;
    uint[][] _column_mat;
    uint[][] _row_mat;


    /**
    各行や各列のどのインデックスに"1"が存在するかの集合を作る
    */
    void generateCompactRep()
    {
        _column_mat = new uint[][](_N);
        _row_mat = new uint[][](_M);

        foreach(i_col; 0 .. _N) {
            foreach(i_row; 0 .. _M) {
                if(_H_mat[i_col][i_row] == 1) {
                    _column_mat[i_col] ~= i_row;
                    _row_mat[i_row] ~= i_col;
                }
            }
        }
    }


    void generateQCLDPC(size_t N, size_t M)(in ref shared int[N][M] baseH)
    {
        _H_mat = new ubyte[][](_N, _M);

        foreach(i_base_row, base_row; baseH)
            foreach(i_base_col, e_base; base_row) if(e_base >= 0) {
                foreach(i; 0 .. _Z)
                    _H_mat[_Z * i_base_col + (i + e_base) % _Z][_Z * i_base_row + i] = 1;
            }

        generateCompactRep();
    }
}


unittest
{
    import std.algorithm;
    import std.random;
    import std.range;
    import std.array;

    BLDPCCode code = new BLDPCCode(648, 324);
    code.loadWiFiLDPC(648, 0);    // 648, 324, 1/2
    
    ubyte[] rndbits = [
        0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,
        0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1,
        1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
        1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0,
        1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1,
        0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0,
        1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
        1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0,
        1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0,
        1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1,
        1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 0
    ];

    auto codeword = code.encode(rndbits);

    assert(rndbits[0 .. 324] == codeword[0 .. 324]);
    assert(code.check_codeword(codeword));
    assert(codeword == [
        0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1,
        0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1,
        1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
        1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
        1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0,
        1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1,
        0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0,
        1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0,
        1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 
        1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0,
        1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1,
        1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0,
        1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1,
        1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1,
        0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
        1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,
        1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0,
        0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1,
        1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1,
        0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1,
        1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1,
        0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0,
        0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0,
        1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1,
        0, 1, 0, 1, 1, 0, 1, 1
    ]);

    import std.math;
    import std.random;
    import simpleldpc.constellation;
    auto mod = new Constellation(1);
    auto syms = mod.modulate(codeword);
    Random rnd;
    rnd.seed(0);
    foreach(ref e; syms) {
        auto x = uniform01(rnd),
             y = uniform01(rnd);

        e += sqrt(-2*log(x))*cos(2*PI*y) * SQRT1_2;
    }
    auto llr = mod.llr_compute(syms, 0.5);
    auto decoded = code.decode(llr, 20);
    assert(decoded == codeword);
}
