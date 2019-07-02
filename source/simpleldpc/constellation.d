module simpleldpc.constellation;

import std.math;
import std.range;
import std.stdio;

class Constellation
{
    this(uint n_bits)
    {
        _n_bits = n_bits;
        _n_syms = 1 << _n_bits;

        // _bit_sym_mapは単にインデックスをbit列に変換しただけ
        _bit_sym_map = new int[][](_n_syms, _n_bits);

        foreach(i_sym; 0 .. _n_syms) {
            uint tmp_sym = i_sym;
            foreach(i_bit; 0 .. _n_bits) {
                if(tmp_sym % 2 == 1)
                    _bit_sym_map[i_sym][i_bit] = 1;
                else
                    _bit_sym_map[i_sym][i_bit] = 0;

                tmp_sym >>= 1;
            }
        }

        // _pointsはグレイコードになるように配置
        switch(_n_bits) {
            case 1:
                _points = [-1, 1];
                break;
            case 2:
                _points = [-3, -1, 3, 1];
                break;
            case 3:
                _points = [-7, -5, -1, -3, 7, 5, 1, 3];
                break;
            default:
                writeln("Constellation not supported");
        }

        double energy = 0;
        foreach(i_sym; 0 .. _n_syms)
            energy += _points[i_sym]^^2;

        energy /= _n_syms;

        foreach(i_sym; 0 .. _n_syms)
            _points[i_sym] /= sqrt(energy);
    }


    double[] modulate(ubyte[] coded_bits)
    {
        auto mod_sym = new double[](coded_bits.length / _n_bits);

        foreach(c_bit; iota(0, coded_bits.length, _n_bits)) {
            // coded_bitsを_n_bitsごとに区切り，シンボルのインデックスへ変換
            uint sym;
            foreach(i_bit; c_bit .. c_bit + _n_bits)
                sym += coded_bits[i_bit] << (i_bit - c_bit);

            mod_sym[c_bit / _n_bits] = _points[sym];
        }

        return mod_sym;
    }


    double[] computeLLR(double[] y, double n_0)
    {
        auto dst = computeP0P1(y, n_0);
        foreach(ref e; dst)
            e = log(e);

        return dst;
    }


    double[] computeP0P1(double[] y, double n_0)
    {
        auto p_0 = new double[](y.length * _n_bits);
        auto p_1 = new double[](y.length * _n_bits);
        auto p0p1 = new double[](y.length * _n_bits);

        p_0[] = 0;
        p_1[] = 0;

        foreach(i_y, e_y; y){
            foreach(i_p, e_p; _points) {
                double p_sym = exp(-(e_y - e_p)^^2 / 2 / n_0);

                foreach(i_bit; 0 .. _n_bits) {
                    if(_bit_sym_map[i_p][i_bit] == 1)
                        p_1[i_y * _n_bits + i_bit] += p_sym;
                    else
                        p_0[i_y * _n_bits + i_bit] += p_sym;
                }
            }
        }

        foreach(i_bit; 0 .. p0p1.length)
            p0p1[i_bit] = p_0[i_bit] / p_1[i_bit];

        return p0p1;
    }


  private:
    double[] _points;
    uint _n_bits;
    uint _n_syms;
    int[][] _bit_sym_map;
}


unittest
{
    auto cs1 = new Constellation(1);
    assert(cs1._bit_sym_map == [[0], [1]]);
    auto syms = cs1.modulate([0, 1]);
    assert(syms.approxEqual([-1, 1]));
    assert(cs1.computeLLR(syms, 1).approxEqual([2, -2]));

    auto cs2 = new Constellation(2);
    assert(cs2._bit_sym_map == [[0, 0], [1, 0], [0, 1], [1, 1]]);
    syms = cs2.modulate([0, 0, 1, 0, 0, 1, 1, 1]);
    assert(syms.approxEqual([-1.34164, -0.447214, 1.34164, 0.447214]));
    assert(cs2.computeLLR(syms, 1).approxEqual([0.163675, 1.98609, -0.649733, 0.649733, 0.163675, -1.98609, -0.649733, -0.649733]));

    auto cs3 = new Constellation(3);
    assert(cs3._bit_sym_map == [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]);
    syms = cs3.modulate([0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1]);
    assert(syms.approxEqual([-1.52753, -1.09109, -0.218218, -0.654654, 1.52753, 1.09109, 0.218218, 0.654654]));
    assert(cs3.computeLLR(syms, 1).approxEqual([
        -0.0389527, 0.319143, 2.14779,
        -0.0798453, -0.0940412, 1.51827,
        -0.0565873, -0.693676, 0.300105,
        -0.0745835, -0.460898, 0.904003,
        -0.0389527, 0.319143, -2.14779,
        -0.0798453, -0.0940412, -1.51827,
        -0.0565873, -0.693676, -0.300105,
        -0.0745835, -0.460898, -0.904003
    ]));
}
