#include "roaring.c"
#include <vector>

class FE_roaring {
    private:
        roaring_bitmap_t **r;
        uint32_t n_bitmaps_;

    public:
        FE_roaring(std::vector<uint32_t>* bitmaps, uint32_t n_keys, uint32_t n_bitmaps) : n_bitmaps_(n_bitmaps) {
            r = (roaring_bitmap_t**)malloc(sizeof(roaring_bitmap_t*) * n_bitmaps);
            for (size_t i = 0; i < n_bitmaps; i++)
            {
                r[i] = roaring_bitmap_create();
                for (size_t j = 0; j < n_keys; j++)
                {
                    roaring_bitmap_add(r[i], bitmaps[i][j]);
                } 
            }
        }

        FE_roaring(uint64_t** bitmaps, uint32_t n_element, uint32_t n_bitmaps) : n_bitmaps_(n_bitmaps) {
            r = (roaring_bitmap_t**)malloc(sizeof(roaring_bitmap_t*) * n_bitmaps);
            for (size_t i = 0; i < n_bitmaps; i++)
            {
                r[i] = roaring_bitmap_create();
                for (size_t j = 0; j < n_element; j++)
                {
                    if ((bitmaps[i][j / 64] & (1ULL << j % 64)) != 0)
                        roaring_bitmap_add(r[i], j);
                } 
            }
        }

        ~FE_roaring() {
            for (size_t i = 0; i < n_bitmaps_; i++)
            {
                roaring_bitmap_free(r[i]);
            }
        }

        size_t read(size_t pos) {
            size_t window_index = 0;
            for (size_t i = 0; i < n_bitmaps_; i++)
            {
                window_index += roaring_bitmap_contains(r[i], pos) ? (1 << i) : 0;
            }
            return window_index;
        }

        void update(size_t pos) {
            for (size_t i = 0; i < n_bitmaps_; i++)
            {
                if (roaring_bitmap_contains(r[i], pos))
                    roaring_bitmap_remove(r[i], pos);
                else {
                    roaring_bitmap_add(r[i], pos);
                    return;
                }
            }
        }

        size_t bitCount() {
            size_t bc = 0;
            for (size_t i = 0; i < n_bitmaps_; i++) {
                bc += (roaring_bitmap_frozen_size_in_bytes(r[i]) << 3);
            }
            return bc + sizeof(*this) * 8 - sizeof(r) * 8;
        }

};