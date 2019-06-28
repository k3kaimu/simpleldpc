module simpleldpc.ldpccode;

import std.exception : enforce;
import std.format : format;
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
    bool check_codeword(in ubyte[] decoded_cw)
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


    uint get_info_length()
    {
        return _K;
    }


    /++
    S. Myung, K. Yang, J. Kim, "Quasi-Cyclic LDPC Codes for Fast Encoding", IEEE Trans. on Inform. Theory, vol. 51, no. 8, August 2005.
    +/
    ubyte[] encode(in ubyte[] info_bits)
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


    // ubyte[] decode(double[] llr, uint max_iter, bool min_sum);

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
}
