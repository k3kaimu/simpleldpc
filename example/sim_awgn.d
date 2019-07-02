/++ dub.json:
{
    "name": "sim_awgn",
    "dependencies": {
        "simpleldpc": { "path": ".." }
    }
}
+/
import simpleldpc.ldpccode;
import simpleldpc.constellation;

import std.algorithm;
import std.array;
import std.datetime;
import std.math;
import std.parallelism;
import std.random;
import std.range;
import std.stdio;



void main()
{    immutable
        N = 648,
        K = 324,
        M = N - K,
        R = K * 1.0 / N,
        Nc = 1;             // 1: QPSK, 2: 16QAM, 3: 64QAM

    BLDPCCode ldpc = new BLDPCCode(N, K);
    ldpc.loadWiFiLDPC(N, 0);    // 648, 324, 1/2

    Constellation mod = new Constellation(Nc);

    immutable
        MaxTotalBlocks = 10_000,        // 1万ブロックまで伝送する
        MaxErrorBlocks  = 100;          // 100ブロック誤った時点で終わり

    auto ebn0_dBList = iota(0, 5.5, 0.25).array;
    double[] berList = new double[ebn0_dBList.length];
    double[] blerList = new double[ebn0_dBList.length];

    foreach(i_ebn0, ebn0_dB; ebn0_dBList.parallel) {
        immutable double ebn0 = 10.0^^(ebn0_dB / 10);
        immutable double snr = ebn0 * K / N * Nc;
        immutable double N0 = 0.5 / snr;
        immutable double noiseAmp = sqrt(N0);

        ubyte[] info = new ubyte[K];
        double[] noise = new double[N / Nc];
        double[] received = new double[N / Nc];

        size_t totalBlocks = 0,
               errorBlocks = 0,
               totalBits = 0,
               errorBits = 0;

        Random rnd;
        rnd.seed(0);
        
        StopWatch sw;
        while(totalBlocks < MaxTotalBlocks && errorBlocks < MaxErrorBlocks) {
            rnd.makeBits(info);

            enum size_t P = 4;
            double[] p0p1;
            ubyte[] coded;
            ubyte[] decoded = new ubyte[N*P];
            
            foreach(i; 0 .. P) {
                coded ~= ldpc.encode(info);
                double[] sym = mod.modulate(coded[i*N .. (i+1)*N]);
                rnd.makeNoise(noise);

                sym[] += noise[] * noiseAmp;

                p0p1 ~= mod.computeP0P1(sym, N0);
                // auto llr = mod.computeLLR(sym, N0);
            }

            sw.start();
            ldpc.decodeP0P1AVX!(float[4])(p0p1, 20, decoded);
            // ubyte[] decoded = ldpc.decodeP0P1(p0p1, 20);
            // ubyte[] decoded = ldpc.decodeLLR(llr, 20);
            sw.stop();

            foreach(i; 0 .. P) {
                if(coded[i*N .. i*N + K] != decoded[i*N .. i*N + K])
                    ++errorBlocks;

                foreach(j; 0 .. K)
                    if(coded[i*N + j] != decoded[i*N + j])
                        ++errorBits;
            }

            totalBits += K * P;
            totalBlocks += P;
        }

        berList[i_ebn0] = errorBits * 1.0 / totalBits;
        blerList[i_ebn0] =  errorBlocks * 1.0 / totalBlocks;
        writefln!"EbN0 = %s [dB], DecodeThroughput = %s [kbps]"(
            ebn0_dB,
            totalBits / (sw.peek.usecs / 1.0E6) / 1E3
        );
    }

    foreach(i, e; ebn0_dBList) {
        writefln!"EbN0 = %s [dB], BER = %s, BLER = %s"(e, berList[i], blerList[i]);
    }
}


void makeNoise(ref Random rnd, double[] dst)
in { assert(dst.length % 2 == 0); }
body {
    foreach(i; 0 .. dst.length / 2) {
        immutable 
            x = uniform01(rnd),
            y = uniform01(rnd);

        immutable r = sqrt(-2 * log(x));
        dst[i*2 + 0] = r * cos(2 * PI * y);
        dst[i*2 + 1] = r * sin(2 * PI * y);
    }
}


void makeBits(ref Random rnd, ubyte[] dst)
{
    foreach(ref e; dst) {
        e = rnd.front % 2;
        rnd.popFront();
    }
}
