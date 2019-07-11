module simpleldpc.ldpccode;

import std.algorithm : min, max;
import std.exception : enforce;
import std.format : format;
import std.math : tanh, atanh;
import std.range : iota;

import simpleldpc.wifildpc;

import core.simd;

version(LDC)
{
    import ldc.attributes;

    pragma(LDC_intrinsic, "llvm.minnum.v4f32")
    float4 _mm_minps_(float4, float4) pure @safe;
    pragma(LDC_intrinsic, "llvm.maxnum.v4f32")
    float4 _mm_maxps_(float4, float4) pure @safe;
    pragma(LDC_intrinsic, "llvm.minnum.v8f32")
    float8 _mm_minps_(float8, float8) pure @safe;
    pragma(LDC_intrinsic, "llvm.maxnum.v8f32")
    float8 _mm_maxps_(float8, float8) pure @safe;
}

version(DigitalMars)
{
    float4 _mm_minps_(float4 x, float4 y) pure @safe { return simd!(XMM.MINPS)(x, y); }
    float4 _mm_maxps_(float4 x, float4 y) pure @safe { return simd!(XMM.MAXPS)(x, y); }
}


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
        return .checkCodeword(_row_mat, decoded_cw);
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
                double product = 1;
                foreach(i_col_index, i_col; row) {
                    // 変数ノードでの計算時に自身のノードのLLRも加えているので
                    // それは引き算して周辺化した値を得る
                    double l1 = updated_llr[i_col] - last_edge_mat[i_row][i_col_index];
                    // tanhの発散を防ぐ
                    l1 = min(l1, 20.0);
                    l1 = max(l1, -20.0);
                    product *= tanh(l1/2);
                }

                foreach(i_col_index1, i_col_1; row) {
                    double l1 = updated_llr[i_col_1] - last_edge_mat[i_row][i_col_index1];
                    l1 = min(l1, 20.0);
                    l1 = max(l1, -20.0);
                    double tmp = product / tanh(l1/2);
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


    ubyte[] decodeP0P1(F)(in F[] input_p0p1, uint max_iter, ubyte[] decoded_cw) const
    {
        static float[][] edge_mat;
        static float[][] last_edge_mat;
        static float[] updated_p0p1;

        if(edge_mat is null) {
            edge_mat = new float[][](_M);
            last_edge_mat = new float[][](_M);
            updated_p0p1 = new float[_N];

            foreach(i; 0 .. _M) {
                edge_mat[i] = new float[](_row_mat[i].length);
                last_edge_mat[i] = new float[](_row_mat[i].length);
            }
        }

        foreach(i; 0 .. _N)
            updated_p0p1[i] = input_p0p1[i];

        foreach(i; 0 .. _M) {
            edge_mat[i][] = 1;
            last_edge_mat[i][] = 1;
        }

        foreach(iter; 0 .. max_iter) 
        {
            foreach(i_row, row; _row_mat) {
                float prob_product = 1;
                foreach(i_col_index, i_col; row) {
                    float p0p1 = updated_p0p1[i_col] / last_edge_mat[i_row][i_col_index];
                    float p1 = 1/(1 + p0p1);
                    prob_product *= (1 - 2*p1);
                }

                foreach(i_col_index1, i_col_1; row) {
                    float p0p1 = updated_p0p1[i_col_1] / last_edge_mat[i_row][i_col_index1];
                    float p1 = 1/(1 + p0p1);
                    auto tmp = prob_product / (1 - 2*p1);

                    tmp = (1 + tmp)/(1 - tmp);
                    tmp = min(tmp, 1E3);
                    tmp = max(tmp, 1E-3);
                    edge_mat[i_row][i_col_index1] = tmp;
                }
            }

            foreach(i; 0 .. _M)
                last_edge_mat[i][] = edge_mat[i][];

            foreach(i; 0 .. _N)
                updated_p0p1[i] = input_p0p1[i];

            foreach(i_row, row; _row_mat) {
                foreach(i_col_index, i_col; row) {
                    updated_p0p1[i_col] *= last_edge_mat[i_row][i_col_index];
                }
            }

            foreach(i; 0 .. _N) {
                if(updated_p0p1[i] > 1)
                    decoded_cw[i] = 0;
                else
                    decoded_cw[i] = 1;
            }

            if(checkCodeword(decoded_cw))
                break;
        }

        return decoded_cw;
    }


    void decodeP0P1SIMD(V, F)(ref Workspace!V ws,in F[] input_p0p1, uint max_iter, ubyte[] decoded_cw) const
    in{
        assert(input_p0p1.length == _N * V.length);
        assert(decoded_cw.length == _N * V.length);
    }
    body {
        sumProductDecodeP0P1SIMD!(V, F)(ws, _row_mat, input_p0p1, max_iter, decoded_cw);
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
    assert(code.checkCodeword(codeword));
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
    auto llr = mod.computeLLR(syms, 0.5);
    auto decoded = code.decodeLLR(llr, 20);
    assert(decoded == codeword);
}


@fastmath
V vecminmax(V, F)(V v, F vmin_, F vmax_) @trusted
{
    V vmin = vmin_,
      vmax = vmax_;

    return _mm_maxps_(_mm_minps_(v, vmax), vmin);
}


struct Workspace(T)
{
    size_t maxRW;

    static if(is(T : float))
        alias E = T;
    else
        alias E = Vector!T;

    E[] edge_mat;
    E[] updated_p0p1;
    E[] input_p0p1_copy;
    E[] prob1m2p;
}


@fastmath
void sumProductDecodeP0P1SIMD(V, F)(ref Workspace!V ws, in uint[][] _row_mat, in F[] input_p0p1, uint max_iter, ubyte[] decoded_cw)
{
    with(ws){
    import core.simd;
    import std.experimental.allocator.mallocator;

    alias alloc = AlignedMallocator.instance;

    static if(is(V : float)) {
        enum size_t P = 1;
        alias VecType = V;
    }else{
        enum size_t P = V.length;
        alias VecType = Vector!V;
    }

    immutable size_t _M = _row_mat.length;
    immutable size_t _N = input_p0p1.length / P;

    void allocate(ref VecType[] v, size_t n) @trusted {
        v = (cast(VecType*)alloc.alignedAllocate(n * VecType.sizeof, VecType.sizeof))[0 .. n];
    }

    if(edge_mat is null) {
        maxRW = 0;
        size_t totElem = 0;
        foreach(row; _row_mat) {
            maxRW = max(maxRW, row.length);
            totElem += row.length;
        }

        allocate(edge_mat, totElem);
        allocate(updated_p0p1, _N);
        allocate(input_p0p1_copy, _N);
        allocate(prob1m2p, maxRW);
    }

    bool[P] success;

    static if(P != 1) {
        foreach(i; 0 .. P) foreach(j; 0 .. _N) {
            input_p0p1_copy[j][i] = input_p0p1[i*_N + j];
            updated_p0p1[j][i] = input_p0p1[i*_N + j];
        }
    } else {
        foreach(i; 0 .. _N) {
            input_p0p1_copy[i] = input_p0p1[i];
            updated_p0p1[i] = input_p0p1[i];
        }
    }

    edge_mat[] = 1;

    foreach(iter; 0 .. max_iter) 
    {
        {
            auto p_last_edge_mat = edge_mat.ptr;
            auto p_edge_mat = edge_mat.ptr;
            foreach(i_row, row; _row_mat) {
                VecType prob_product = 1;
                foreach(i_col_index1, i_col_1; row) {
                    VecType q1 = updated_p0p1[i_col_1];
                    VecType q2 = *p_last_edge_mat;
                    ++p_last_edge_mat;

                    VecType p1m2p = (q1 - q2) / (q1 + q2);
                    prob1m2p[i_col_index1] = p1m2p;
                    prob_product *= p1m2p;
                }

                foreach(i_col_index1, i_col_1; row) {
                    VecType p1m2p = prob1m2p[i_col_index1];
                    VecType tmp = (p1m2p + prob_product) / (p1m2p - prob_product);
                    *p_edge_mat = vecminmax(tmp, 0.5^^15, 2.0^^15);
                    ++p_edge_mat;
                }
            }
        }

        updated_p0p1[] = input_p0p1_copy[];

        {
            auto p_edge_mat = edge_mat.ptr;
            foreach(i_row, row; _row_mat) {
                foreach(i_col_index, i_col; row) {
                    updated_p0p1[i_col] *= *p_edge_mat;
                    ++p_edge_mat;
                }
            }
        }

        foreach(i; 0 .. _N) {
            updated_p0p1[i] = vecminmax(updated_p0p1[i], 0, 2.0^^25);
        }

        foreach(i; 0 .. P){
            if(success[i]) continue;

            foreach(j; 0 .. _N) {
                if(updated_p0p1[j][i] > 1)
                    decoded_cw[i*_N + j] = 0;
                else
                    decoded_cw[i*_N + j] = 1;
            }

            if(checkCodeword(_row_mat, decoded_cw[i*_N .. (i+1)*_N]))
                success[i] = true;
        }

        bool checkAllSuccess = true;
        foreach(i; 0 .. P)
            checkAllSuccess = checkAllSuccess && success[i];

        if(checkAllSuccess)
            break;
    }
    }
}


bool checkCodeword(in uint[][] _row_mat, in ubyte[] decoded_cw) @safe
{
    immutable _M = _row_mat.length;

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
