/* Fixed switching data from Simoens--Van Overberghe,
 * "An algorithm to find cospectral mates" (2026):
 * https://github.com/robinsimoens/implementing-switching-methods
 *
 * The irreducible switching-subgraph catalogues were adapted with the
 * authors' permission from their SageMath reference implementation and
 * re-expressed here as compact edge masks under this project's GPL-3.0
 * license. */

#include "switching_methods.h"

#include <stdint.h>

static const int8_t gm4_numerator[] = {
        -1, 1, 1, 1, 1, -1, 1, 1, 1, 1,
        -1, 1, 1, 1, 1, -1,
};

static const uint64_t gm4_subgraphs[] = {
        UINT64_C(0x0), UINT64_C(0x21), UINT64_C(0x2d), UINT64_C(0x3f),
};

static const int8_t gm6_numerator[] = {
        -2, 1, 1, 1, 1, 1, 1, -2, 1, 1,
        1, 1, 1, 1, -2, 1, 1, 1, 1, 1,
        1, -2, 1, 1, 1, 1, 1, 1, -2, 1,
        1, 1, 1, 1, 1, -2,
};

static const uint64_t gm6_subgraphs[] = {
        UINT64_C(0x0), UINT64_C(0x884), UINT64_C(0xd4c), UINT64_C(0xfdc),
        UINT64_C(0x254a), UINT64_C(0x2dce), UINT64_C(0x3dfe), UINT64_C(0x7fff),
};

static const int8_t wqh6_numerator[] = {
        2, -1, -1, 1, 1, 1, -1, 2, -1, 1,
        1, 1, -1, -1, 2, 1, 1, 1, 1, 1,
        1, 2, -1, -1, 1, 1, 1, -1, 2, -1,
        1, 1, 1, -1, -1, 2,
};

static const uint64_t wqh6_subgraphs[] = {
        UINT64_C(0x0), UINT64_C(0x2212), UINT64_C(0x6822), UINT64_C(0x6732),
        UINT64_C(0x884), UINT64_C(0x6b2a), UINT64_C(0xd4c), UINT64_C(0x4a9a),
        UINT64_C(0xfdc), UINT64_C(0x7023), UINT64_C(0x4c5a), UINT64_C(0x6fb6),
        UINT64_C(0x732b), UINT64_C(0x4f9e), UINT64_C(0x777b), UINT64_C(0x7fff),
};

static const int8_t ah6_numerator[] = {
        1, 1, 0, 0, 1, -1, 1, 1, 0, 0,
        -1, 1, 1, -1, 1, 1, 0, 0, -1, 1,
        1, 1, 0, 0, 0, 0, 1, -1, 1, 1,
        0, 0, -1, 1, 1, 1,
};

static const uint64_t ah6_subgraphs[] = {
        UINT64_C(0x3116), UINT64_C(0x7371),
};

static const int8_t is6_numerator[] = {
        2, 3, 3, -1, -1, -1, 3, 2, -3, 1,
        1, 1, 3, -3, 2, 1, 1, 1, -1, 1,
        1, -2, 3, 3, -1, 1, 1, 3, -2, 3,
        -1, 1, 1, 3, 3, -2,
};

static const uint64_t is6_subgraphs[] = {
        UINT64_C(0x0), UINT64_C(0xe3d), UINT64_C(0x66d3), UINT64_C(0x192c),
        UINT64_C(0x71c2), UINT64_C(0x7fff),
};

static const int8_t fano_numerator[] = {
        -1, 1, 1, 0, 1, 0, 0, 0, -1, 1,
        1, 0, 1, 0, 0, 0, -1, 1, 1, 0,
        1, 1, 0, 0, -1, 1, 1, 0, 0, 1,
        0, 0, -1, 1, 1, 1, 0, 1, 0, 0,
        -1, 1, 1, 1, 0, 1, 0, 0, -1,
};

static const uint64_t fano_subgraphs[] = {
        UINT64_C(0x107b1a), UINT64_C(0x2630c), UINT64_C(0x1eacd), UINT64_C(0xb779e),
};

static const int8_t gm8_numerator[] = {
        -3, 1, 1, 1, 1, 1, 1, 1, 1, -3,
        1, 1, 1, 1, 1, 1, 1, 1, -3, 1,
        1, 1, 1, 1, 1, 1, 1, -3, 1, 1,
        1, 1, 1, 1, 1, 1, -3, 1, 1, 1,
        1, 1, 1, 1, 1, -3, 1, 1, 1, 1,
        1, 1, 1, 1, -3, 1, 1, 1, 1, 1,
        1, 1, 1, -3,
};

static const uint64_t gm8_subgraphs[] = {
        UINT64_C(0x0), UINT64_C(0x210408), UINT64_C(0x330618), UINT64_C(0x328a18),
        UINT64_C(0x3b5638), UINT64_C(0x3fde78), UINT64_C(0x10b0a14), UINT64_C(0x1338e1c),
        UINT64_C(0x734e54), UINT64_C(0x72ce64), UINT64_C(0x6b5a34), UINT64_C(0x13bde3c),
        UINT64_C(0x517dd3c), UINT64_C(0x497dd5c), UINT64_C(0x59bff7c), UINT64_C(0x4a9552a),
        UINT64_C(0x4bb573a), UINT64_C(0x4badb3a), UINT64_C(0x4bfdf7a), UINT64_C(0x5bbdf3e),
        UINT64_C(0x7bfdffe), UINT64_C(0xfffffff),
};

static const int8_t gm44_numerator[] = {
        -1, 1, 1, 1, 0, 0, 0, 0, 1, -1,
        1, 1, 0, 0, 0, 0, 1, 1, -1, 1,
        0, 0, 0, 0, 1, 1, 1, -1, 0, 0,
        0, 0, 0, 0, 0, 0, -1, 1, 1, 1,
        0, 0, 0, 0, 1, -1, 1, 1, 0, 0,
        0, 0, 1, 1, -1, 1, 0, 0, 0, 0,
        1, 1, 1, -1,
};

static const uint64_t gm44_subgraphs[] = {
        UINT64_C(0x210408), UINT64_C(0x85122), UINT64_C(0x211c68), UINT64_C(0x851a6),
        UINT64_C(0x4a11c68), UINT64_C(0x7a11c68), UINT64_C(0x36e241), UINT64_C(0x10a3c7),
        UINT64_C(0x36e343), UINT64_C(0x36e3c7), UINT64_C(0x4885122), UINT64_C(0x3085122),
        UINT64_C(0x48851a6), UINT64_C(0x4bb5122), UINT64_C(0x33b5122), UINT64_C(0x4bb51a6),
        UINT64_C(0x3b5638), UINT64_C(0x78851a6), UINT64_C(0x4bb5638), UINT64_C(0x7bb5638),
        UINT64_C(0x7bb51a6), UINT64_C(0xffb5638), UINT64_C(0x7860894), UINT64_C(0x7b608f4),
        UINT64_C(0x4b608f4), UINT64_C(0x4ae1d12), UINT64_C(0xfc60912), UINT64_C(0xcf60972),
        UINT64_C(0xcf51172), UINT64_C(0xff60972), UINT64_C(0x7860996), UINT64_C(0x7b609f6),
        UINT64_C(0xcf61d22), UINT64_C(0xff61d22), UINT64_C(0x7bb71a7), UINT64_C(0xcf609f6),
        UINT64_C(0xcf511f6), UINT64_C(0xff609f6), UINT64_C(0xfc60996), UINT64_C(0x4bb573a),
        UINT64_C(0xfc871a7), UINT64_C(0x4adddee), UINT64_C(0x4badb3a), UINT64_C(0x7badb3a),
        UINT64_C(0xffadb3a), UINT64_C(0xffb71a7), UINT64_C(0x7adddee), UINT64_C(0x7af5dbe),
        UINT64_C(0x79efbf7), UINT64_C(0xfdefbf7),
};

static const int8_t wqh8_numerator[] = {
        3, -1, -1, -1, 1, 1, 1, 1, -1, 3,
        -1, -1, 1, 1, 1, 1, -1, -1, 3, -1,
        1, 1, 1, 1, -1, -1, -1, 3, 1, 1,
        1, 1, 1, 1, 1, 1, 3, -1, -1, -1,
        1, 1, 1, 1, -1, 3, -1, -1, 1, 1,
        1, 1, -1, -1, 3, -1, 1, 1, 1, 1,
        -1, -1, -1, 3,
};

static const uint64_t wqh8_subgraphs[] = {
        UINT64_C(0x0), UINT64_C(0x1040044), UINT64_C(0x50c1144), UINT64_C(0xd1e3144),
        UINT64_C(0xd200146), UINT64_C(0x9241124), UINT64_C(0xd2e3124), UINT64_C(0x48c1124),
        UINT64_C(0xc9e3124), UINT64_C(0x210408), UINT64_C(0x4800102), UINT64_C(0x90610a4),
        UINT64_C(0x5840146), UINT64_C(0x10c5064), UINT64_C(0xc1e1954), UINT64_C(0xd962146),
        UINT64_C(0x8250854), UINT64_C(0xd2a6066), UINT64_C(0xc2e1934), UINT64_C(0xda62126),
        UINT64_C(0x8251910), UINT64_C(0xd2a70a4), UINT64_C(0x52cf064), UINT64_C(0xc2d1954),
        UINT64_C(0x78e2126), UINT64_C(0x330618), UINT64_C(0x7800186), UINT64_C(0x4885122),
        UINT64_C(0x58c5166), UINT64_C(0x9966066), UINT64_C(0xd9e7166), UINT64_C(0x328a18),
        UINT64_C(0x3085122), UINT64_C(0xd8611a6), UINT64_C(0x90e5962), UINT64_C(0xb1a7122),
        UINT64_C(0xf922186), UINT64_C(0x81ef860), UINT64_C(0xb1a6942), UINT64_C(0x9a44966),
        UINT64_C(0xdae6966), UINT64_C(0xb284966), UINT64_C(0xd2ef962), UINT64_C(0xdad7166),
        UINT64_C(0xcb0e066), UINT64_C(0xe3f3534), UINT64_C(0xcaa6922), UINT64_C(0xd8e59e6),
        UINT64_C(0x3b5638), UINT64_C(0x78851a6), UINT64_C(0x48cd962), UINT64_C(0xf9a71a6),
        UINT64_C(0xc9ef962), UINT64_C(0xfaa69a6), UINT64_C(0x3fde78), UINT64_C(0x78cd9e6),
        UINT64_C(0xf9ef9e6), UINT64_C(0x8308264), UINT64_C(0x52c5954), UINT64_C(0x52d154c),
        UINT64_C(0xd3f354c), UINT64_C(0x50e5952), UINT64_C(0x9165952), UINT64_C(0x49e152c),
        UINT64_C(0x82f0c6c), UINT64_C(0xb250566), UINT64_C(0xcbe7934), UINT64_C(0x5a56056),
        UINT64_C(0xdac7965), UINT64_C(0xca50566), UINT64_C(0xcbf352c), UINT64_C(0x42f1d18),
        UINT64_C(0xc871992), UINT64_C(0xc3f1d5c), UINT64_C(0xdb66956), UINT64_C(0x7860996),
        UINT64_C(0x48e5932), UINT64_C(0x79e7136), UINT64_C(0x79e6956), UINT64_C(0xb2d5976),
        UINT64_C(0xfa86967), UINT64_C(0xe2f0d76), UINT64_C(0xc3fbd48), UINT64_C(0xdb67994),
        UINT64_C(0xd9659d6), UINT64_C(0x48d5952), UINT64_C(0x48d58d4), UINT64_C(0x48cf861),
        UINT64_C(0xcaf78b4), UINT64_C(0xb2f78b4), UINT64_C(0x8368c6c), UINT64_C(0xe2fa876),
        UINT64_C(0x78e59b6), UINT64_C(0x78cf963), UINT64_C(0xfaf79b6), UINT64_C(0xcb4d976),
        UINT64_C(0xf3ef976), UINT64_C(0x83f9a7c), UINT64_C(0xcb64d76), UINT64_C(0xfb0a367),
        UINT64_C(0xb364d76), UINT64_C(0x4bd9b5c), UINT64_C(0xb4e59b6), UINT64_C(0xf8fbbe3),
        UINT64_C(0xe3ffd72), UINT64_C(0x4a9552a), UINT64_C(0xfc02187), UINT64_C(0xc3b5d5a),
        UINT64_C(0xcb68b76), UINT64_C(0x9b5556e), UINT64_C(0xe3ba276), UINT64_C(0xf3fb376),
        UINT64_C(0x4bb573a), UINT64_C(0x7a955ae), UINT64_C(0xfc871a7), UINT64_C(0x4badb3a),
        UINT64_C(0x79693b6), UINT64_C(0xdaf5dee), UINT64_C(0xfbb75ae), UINT64_C(0xcbffd6a),
        UINT64_C(0x4bfdf7a), UINT64_C(0x7adddee), UINT64_C(0xfccf9e7), UINT64_C(0xfbffdee),
        UINT64_C(0x7bcfb75), UINT64_C(0x7af5dbe), UINT64_C(0xfce79b7), UINT64_C(0x7bfdffe),
        UINT64_C(0xfdefbf7), UINT64_C(0xfffffff),
};

static const int8_t is8_level3_numerator[] = {
        2, 1, 0, 0, 1, -1, -1, 1, 1, 2,
        0, 0, -1, 1, 1, -1, 0, 0, 2, 1,
        -1, 1, -1, 1, 0, 0, 1, 2, 1, -1,
        1, -1, 1, -1, -1, 1, 1, 2, 0, 0,
        -1, 1, 1, -1, 2, 1, 0, 0, -1, 1,
        -1, 1, 0, 0, 1, 2, 1, -1, 1, -1,
        0, 0, 2, 1,
};

static const uint64_t is8_level3_subgraphs[] = {
        UINT64_C(0x449021), UINT64_C(0xf920211), UINT64_C(0x8121a71), UINT64_C(0xf921a71),
        UINT64_C(0x1a5428), UINT64_C(0x79a5428), UINT64_C(0xf8869c6), UINT64_C(0xfbb69c6),
        UINT64_C(0x81ec211), UINT64_C(0x81ec397), UINT64_C(0x865aa51), UINT64_C(0xf9eda71),
        UINT64_C(0x865abd7), UINT64_C(0x6dfc68), UINT64_C(0x7a58bd6), UINT64_C(0x7f791a7),
        UINT64_C(0xfe5abd7), UINT64_C(0xf9edbf7),
};

static const int8_t is8_level5_numerator[] = {
        3, 1, 2, -1, -2, 1, 2, -1, -1, 3,
        1, 2, -1, -2, 1, 2, 2, -1, 3, 1,
        2, -1, -2, 1, 1, 2, -1, 3, 1, 2,
        -1, -2, -2, 1, 2, -1, 3, 1, 2, -1,
        -1, -2, 1, 2, -1, 3, 1, 2, 2, -1,
        -2, 1, 2, -1, 3, 1, 1, 2, -1, -2,
        1, 2, -1, 3,
};

static const uint64_t is8_level5_subgraphs[] = {
        UINT64_C(0x0), UINT64_C(0x8000841), UINT64_C(0x14028), UINT64_C(0x1401450),
        UINT64_C(0x2a9450), UINT64_C(0xb428050), UINT64_C(0xb6a9450), UINT64_C(0x3044822),
        UINT64_C(0x950a0a), UINT64_C(0xbd50a0a), UINT64_C(0xb40a014), UINT64_C(0xb68b414),
        UINT64_C(0xf2a8641), UINT64_C(0x18058e2), UINT64_C(0xe510b11), UINT64_C(0xd4a1751),
        UINT64_C(0xcda734c), UINT64_C(0x210408), UINT64_C(0x85020), UINT64_C(0x11010b0),
        UINT64_C(0x238458), UINT64_C(0xad070), UINT64_C(0x7a08442), UINT64_C(0xb638458),
        UINT64_C(0xb4ad070), UINT64_C(0x5389421), UINT64_C(0xf108071), UINT64_C(0xd12231c),
        UINT64_C(0x91a321c), UINT64_C(0xf389471), UINT64_C(0x14ad020), UINT64_C(0xa3038d0),
        UINT64_C(0xb1388a), UINT64_C(0x34a5654), UINT64_C(0xba3c8ea), UINT64_C(0xed91b13),
        UINT64_C(0x844a841), UINT64_C(0x295428), UINT64_C(0x4885122), UINT64_C(0x2bd478),
        UINT64_C(0x17ca78), UINT64_C(0xb42a0d5), UINT64_C(0xd5a9131), UINT64_C(0x7a8d462),
        UINT64_C(0x93e9478), UINT64_C(0xb4e8655), UINT64_C(0xa381c55), UINT64_C(0xaa5e8ea),
        UINT64_C(0xb6bd478), UINT64_C(0xb57ca78), UINT64_C(0x3520085), UINT64_C(0xa350e58),
        UINT64_C(0xa1c5a70), UINT64_C(0xed04b44), UINT64_C(0x9d4384a), UINT64_C(0x4b18648),
        UINT64_C(0x498d260), UINT64_C(0x8b38ac8), UINT64_C(0xf28a605), UINT64_C(0x37a1485),
        UINT64_C(0x41ab251), UINT64_C(0x27e14c1), UINT64_C(0xd483715), UINT64_C(0xf5ab251),
        UINT64_C(0x903ca78), UINT64_C(0x8db9258), UINT64_C(0xa295478), UINT64_C(0xf62a5d5),
        UINT64_C(0xf62654c), UINT64_C(0xff18648), UINT64_C(0xfd8d260), UINT64_C(0xec93393),
        UINT64_C(0x92bde78), UINT64_C(0x378e58), UINT64_C(0x1eda70), UINT64_C(0x4b9d668),
        UINT64_C(0xb778e58), UINT64_C(0xb5eda70), UINT64_C(0xff9d668), UINT64_C(0x3fde78),
        UINT64_C(0xb56aad5), UINT64_C(0xb7fde78), UINT64_C(0xb7ebed5), UINT64_C(0x2204c44),
        UINT64_C(0x2185034), UINT64_C(0x89502a), UINT64_C(0x5c01552), UINT64_C(0xb32208c),
        UINT64_C(0x9d5a2a), UINT64_C(0x3268c52), UINT64_C(0x3254c2a), UINT64_C(0xe285570),
        UINT64_C(0xf42c170), UINT64_C(0xf6ad570), UINT64_C(0x3543812), UINT64_C(0x355783a),
        UINT64_C(0xb185a34), UINT64_C(0xb344c6c), UINT64_C(0xcd481e1), UINT64_C(0xf68f534),
        UINT64_C(0x353c03c), UINT64_C(0x1a15cea), UINT64_C(0xd826b4e), UINT64_C(0x8dc90e1),
        UINT64_C(0xcb8556c), UINT64_C(0x98a7a4e), UINT64_C(0x77ad534), UINT64_C(0x4d1311a),
        UINT64_C(0x9a91e0a), UINT64_C(0xd885b22), UINT64_C(0x3305cf4), UINT64_C(0x64547f1),
        UINT64_C(0x24d56f1), UINT64_C(0x3eadc46), UINT64_C(0x1338e1c), UINT64_C(0xd7b9539),
        UINT64_C(0x37bd43c), UINT64_C(0xf2bc669), UINT64_C(0xdc8d774), UINT64_C(0x9c9d67c),
        UINT64_C(0x79172fa), UINT64_C(0xbfd1e0a), UINT64_C(0xfdc5b22), UINT64_C(0x5f1659),
        UINT64_C(0x7f99e03), UINT64_C(0x11ada34), UINT64_C(0x9e9a52), UINT64_C(0xc5b3319),
        UINT64_C(0xbde9a52), UINT64_C(0xdfad56c), UINT64_C(0xfd66b4e), UINT64_C(0xbf3b8da),
        UINT64_C(0xbcbd07a), UINT64_C(0x63fc669), UINT64_C(0x7bcde62), UINT64_C(0xd5eb3b1),
        UINT64_C(0xf3cb6f1), UINT64_C(0x13bde3c), UINT64_C(0xbf9e5a), UINT64_C(0xbff9e5a),
        UINT64_C(0xfba71ae), UINT64_C(0xbdfda7a), UINT64_C(0x7e2cd46), UINT64_C(0x9ce91b2),
        UINT64_C(0x96785b8), UINT64_C(0xa8d72fa), UINT64_C(0x29f52fa), UINT64_C(0xf572b1b),
        UINT64_C(0x5b1be5a), UINT64_C(0xa4eea45), UINT64_C(0x7dad136), UINT64_C(0x45f5779),
        UINT64_C(0xbeeddf6), UINT64_C(0x4a9552a), UINT64_C(0x4968b52), UINT64_C(0x3266c4e),
        UINT64_C(0x43bb659), UINT64_C(0xf53a259), UINT64_C(0xf7bb659), UINT64_C(0x9c359fa),
        UINT64_C(0xfb44d6e), UINT64_C(0xea9557a), UINT64_C(0x39b52be), UINT64_C(0xf52c975),
        UINT64_C(0xd83cb7a), UINT64_C(0xfc3c17a), UINT64_C(0xfebd57a), UINT64_C(0xfd7cb7a),
        UINT64_C(0xff3f9fa), UINT64_C(0x4b78f5a), UINT64_C(0x49edb72), UINT64_C(0xf5bf279),
        UINT64_C(0x5b6cd6e), UINT64_C(0x59adb36), UINT64_C(0xb6f5d3a), UINT64_C(0xdabdf7a),
        UINT64_C(0xff78f5a), UINT64_C(0xfdedb72), UINT64_C(0x4bfdf7a), UINT64_C(0xb77aedd),
        UINT64_C(0xb5efaf5), UINT64_C(0xfffdf7a), UINT64_C(0xf7efff5), UINT64_C(0xace7bd6),
        UINT64_C(0x9f1bedc), UINT64_C(0x5bbdf3e), UINT64_C(0x73ffe7b), UINT64_C(0x3f5fcfe),
        UINT64_C(0xb7ffefd), UINT64_C(0xfdefbf7), UINT64_C(0xfffffff),
};

static const int8_t ah10_numerator[] = {
        1, 1, 0, 0, 0, 0, 0, 0, 1, -1,
        1, 1, 0, 0, 0, 0, 0, 0, -1, 1,
        1, -1, 1, 1, 0, 0, 0, 0, 0, 0,
        -1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 1, -1, 1, 1, 0, 0, 0, 0,
        0, 0, -1, 1, 1, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 1, -1, 1, 1, 0, 0,
        0, 0, 0, 0, -1, 1, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 1, -1, 1, 1,
        0, 0, 0, 0, 0, 0, -1, 1, 1, 1,
};

static const uint64_t ah10_subgraphs[] = {
        UINT64_C(0x602c805840), UINT64_C(0x1f9fd37fa7bf), UINT64_C(0x602c805e46), UINT64_C(0x1f9fd37fa1b9),
        UINT64_C(0x602f8c5e46), UINT64_C(0x79ac805e46), UINT64_C(0x1f9fd073a1b9), UINT64_C(0x1f86537fa1b9),
        UINT64_C(0x79af8c5e46), UINT64_C(0xf602f8c5e46), UINT64_C(0x1f865073a1b9), UINT64_C(0x109fd073a1b9),
        UINT64_C(0x79af8ddfc6), UINT64_C(0x1f8650722039), UINT64_C(0xf79af8ddfc6), UINT64_C(0x108650722039),
        UINT64_C(0x970ae855cc2), UINT64_C(0x168f517aa33d), UINT64_C(0x970ae855ac4), UINT64_C(0x970ae8544da),
        UINT64_C(0x168f517aa53b), UINT64_C(0x168f517abb25), UINT64_C(0x970ae853aa4), UINT64_C(0x9692e855ac4),
        UINT64_C(0x9709e455ac4), UINT64_C(0x970ae8524ba), UINT64_C(0x168f517ac55b), UINT64_C(0x1696d17aa53b),
        UINT64_C(0x168f61baa53b), UINT64_C(0x168f517adb45), UINT64_C(0x9169e455ac4), UINT64_C(0x16e961baa53b),
        UINT64_C(0xc5c3b217710), UINT64_C(0x13a3c4de88ef), UINT64_C(0xc5c3b217116), UINT64_C(0x13a3c4de8ee9),
        UINT64_C(0xc5c3b21690e), UINT64_C(0xc5c382d7116), UINT64_C(0x13a3c4de96f1), UINT64_C(0x13a3c7d28ee9),
        UINT64_C(0xc5c382d690e), UINT64_C(0xc5c3711690e), UINT64_C(0x13a3c7d296f1), UINT64_C(0x13a3c8ee96f1),
        UINT64_C(0xc5c341d690e), UINT64_C(0x13a3cbe296f1), UINT64_C(0xc5c341ce88e), UINT64_C(0x13a3cbe31771),
};

static const struct switching_method_definition methods[] = {
        {
                GRAPHSWITCHING_METHOD_GM, "gm4",
                4, 2, gm4_numerator,
                gm4_subgraphs,
                sizeof(gm4_subgraphs) / sizeof(gm4_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_GM6, "gm6",
                6, 3, gm6_numerator,
                gm6_subgraphs,
                sizeof(gm6_subgraphs) / sizeof(gm6_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_WQH6, "wqh6",
                6, 3, wqh6_numerator,
                wqh6_subgraphs,
                sizeof(wqh6_subgraphs) / sizeof(wqh6_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_AH6, "ah6",
                6, 2, ah6_numerator,
                ah6_subgraphs,
                sizeof(ah6_subgraphs) / sizeof(ah6_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_IS6, "is6",
                6, 5, is6_numerator,
                is6_subgraphs,
                sizeof(is6_subgraphs) / sizeof(is6_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_FANO, "fano",
                7, 2, fano_numerator,
                fano_subgraphs,
                sizeof(fano_subgraphs) / sizeof(fano_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_GM8, "gm8",
                8, 4, gm8_numerator,
                gm8_subgraphs,
                sizeof(gm8_subgraphs) / sizeof(gm8_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_GM44, "gm44",
                8, 2, gm44_numerator,
                gm44_subgraphs,
                sizeof(gm44_subgraphs) / sizeof(gm44_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_WQH8, "wqh8",
                8, 4, wqh8_numerator,
                wqh8_subgraphs,
                sizeof(wqh8_subgraphs) / sizeof(wqh8_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_IS8_LEVEL3, "is8-level3",
                8, 3, is8_level3_numerator,
                is8_level3_subgraphs,
                sizeof(is8_level3_subgraphs) / sizeof(is8_level3_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_IS8_LEVEL5, "is8-level5",
                8, 5, is8_level5_numerator,
                is8_level5_subgraphs,
                sizeof(is8_level5_subgraphs) / sizeof(is8_level5_subgraphs[0])
        },
        {
                GRAPHSWITCHING_METHOD_AH10, "ah10",
                10, 2, ah10_numerator,
                ah10_subgraphs,
                sizeof(ah10_subgraphs) / sizeof(ah10_subgraphs[0])
        },
};

const struct switching_method_definition *
graphswitching_method_definition(enum graphswitching_method method)
{
        size_t index;

        for (index = 0; index < sizeof(methods) / sizeof(methods[0]);
             ++index) {
                if (methods[index].method == method) {
                        return &methods[index];
                }
        }
        return NULL;
}
