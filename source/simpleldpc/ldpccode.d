module simpleldpc.ldpccode;


class LDPCCode
{
    this(uint block_length, uint info_length)
    {
        _N = block_length;
        _K = info_length;
        _Z = 0;

        _M = _N - _K;
        _H_mat = new ubyte[][](_N, _M);
    }


    // void generate_gallagher_ldpc();

    // bool check_codeword(ubyte[]);

    // void load_wifi_ldpc(uint block_length, uint rate_index);

    // uint get_info_length()
    // {
    //     return _K;
    // }

    // ubyte[] encode(ubyte[] info_bits);

    // ubyte[] decode(double[] llr, uint max_iter, bool min_sum);

  private:
    ubyte[][] _H_mat;
    uint _N;
    uint _K;
    uint _M;
    uint _Z;
    uint[][] _column_mat;
    uint[][] _row_mat;


    // void generate_compact_rep();

    // void lifted_ldpc(int[][] baseH);
}