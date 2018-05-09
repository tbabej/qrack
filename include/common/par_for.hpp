#pragma once

#include <atomic>
#include <future>

#if 1
#define PAR_FOR(BEGIN, END) par_for(BEGIN, END, [&](const bitCapInt idx, const int cpu)
#define PAR_FOR_END );

#define PAR_FOR_SKIP(BEGIN, END, SKIPMASK, MASKWIDTH) par_for_skip(BEGIN, END, SKIPMASK, MASKWIDTH, [&](const bitCapInt idx, const int cpu)
#define PAR_FOR_SKIP_END );

#define PAR_FOR_MASK(BEGIN, END, MASKARRAY, MASKLEN) par_for_mask(BEGIN, END, MASKARRAY, MASKLEN, [&](const bitCapInt idx, const int cpu)
#define PAR_FOR_MASK_END );
#else
#define INC_EQ(I) (I)
#define PAR_FOR_(BEGIN, END, INCFN)\
{\
    std::atomic<bitCapInt> _tidx;\
    _tidx = BEGIN;\
    bitCapInt _pStride = 8;\
\
    std::vector<std::future<void>> _futures(numCores);\
    for (int cpu = 0; cpu < numCores; cpu++) {\
        _futures[cpu] = std::async(std::launch::async, [&, cpu]() {\
            bitCapInt _j, idx = 0;\
            bitCapInt strideEnd = END / _pStride;\
            for (bitCapInt _i = _tidx++; _i < strideEnd; _i = _tidx++) {\
                if (idx >= END) {\
                    break;\
                }\
                for (_j = 0; _j < _pStride; _j++) {\
                    idx = INCFN(_i * _pStride + _j);\
                    /* Easiest to clamp on END. */\
                    if (idx >= END) {\
                        break;\
                    }

#define PAR_FOR_END \
                }\
            }\
        });\
    }\
\
    for (int cpu = 0; cpu < numCores; cpu++) {\
        _futures[cpu].get();\
    }\
}
#define PAR_FOR(BEGIN, END) PAR_FOR_(BEGIN, END, INC_EQ)

#define INC_SKIP(I) ((((I) << _maskWidth) & _highMask) | ((I) & _lowMask))
#define PAR_FOR_SKIP(BEGIN, END, SKIPMASK, MASKWIDTH) \
{\
    /*\
     * Add _maskWidth bits by shifting the incrementor up that number of\
     * bits, filling with 0's.\
     *\
     * For example, if the SKIPMASK is 0x8, then the _lowMask will be 0x7\
     * and the high mask will be ~(0x7 + 0x8) ==> ~0xf, shifted by the\
     * number of extra bits to add.\
     */\
    bitCapInt _maskWidth = MASKWIDTH;\
    bitCapInt _lowMask = (SKIPMASK) - 1;\
    bitCapInt _highMask = (~(_lowMask + (SKIPMASK))) << (_maskWidth - 1);\
\
    PAR_FOR(BEGIN, END, INC_SKIP)

#define PAR_FOR_SKIP_END \
    PAR_FOR_END     \
}

#define INC_MASK(M) ({ \
        bitCapInt _mask = M;\
        for (int _m = 0; _m < _maskLen; _m++) {\
            _mask = ((_mask << 1) & _masks[_m][1]) | (_mask & _masks[_m][0]);\
        }\
        _mask;\
    })
#define PAR_FOR_MASK(BEGIN, END, MASKARRAY, MASKLENGTH) \
{\
    bitLenInt _maskLen = MASKLENGTH;\
    for (int _i = 1; _i < _maskLen; _i++) {\
        if (MASKARRAY[_i] < MASKARRAY[_i - 1]) {\
            throw std::invalid_argument("Masks must be ordered by size");\
        }\
    }\
\
    /* Pre-calculate the masks to simplify the increment function later. */\
    bitCapInt _masks[_maskLen][2];\
\
    for (int _i = 0; _i < _maskLen; _i++) {\
        _masks[_i][0] = MASKARRAY[_i] - 1;  /* low mask */ \
        _masks[_i][1] = (~(_masks[_i][0] + MASKARRAY[_i])); /* high mask */ \
    }\
\
    PAR_FOR(BEGIN, END, INC_MASK)

#define PAR_FOR_MASK_END \
    PAR_FOR_END     \
}
#endif
