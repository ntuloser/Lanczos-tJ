/*
 * This is an optimized version of the following C++ program:
 *
 *   http://keithlea.com/javabench/src/cpp/hash.cpp
 *
 * Keith in his benchmark (http://keithlea.com/javabench/data) showed that the
 * Java implementation is twice as fast as the C++ version. In fact, this is
 * only because the C++ implementation is substandard. Most importantly, Keith
 * is using "sprintf()" to convert an integer to a string, which is known to be
 * extremely inefficient.
 */
#include<iostream>
#include <stdio.h>
#include "khash.h"

KHASH_MAP_INIT_INT64(hMap, double)
int main() {
    int ret;
    int is_missing;
    khiter_t k;
    khash_t(hMap) *h = kh_init(hMap);
    k = kh_put(hMap, h, 5000000000, &ret);
    kh_value(h, k) = 10.232323;
    k = kh_get(hMap, h, 5000000000);
    is_missing = (k == kh_end(h));
    std::cout<<" end ?"<<is_missing<<std::endl;
    k = kh_get(hMap, h, 5);
    kh_del(hMap, h, k);
    for (k = kh_begin(h); k != kh_end(h); ++k)
        if (kh_exist(h, k)) kh_value(h, k) = 1;
    kh_destroy(hMap, h);
    return 0;
}

/*
int main(int argc, char *argv[])
{
    int i, l, n = 1000;
    double ret;
    khash_t(long long) *h, *h2;
    khint_t k;
    h = kh_init(str);
    h2 = kh_init(str);
    if (argc > 1) n = atoi(argv[1]);
    for (i = 0; i < 10000; ++i) {
        char buf[32];
        strcpy(buf, "foo_");
        int2str(i, 10, buf+4);
        k = kh_put(long long, h, strdup(buf), &ret);
        kh_val(h, k) = i;
    }
    for (i = 0; i < n; ++i) {
        for (k = kh_begin(h); k != kh_end(h); ++k) {
            if (kh_exist(h, k)) {
                khint_t k2 = kh_put(str, h2, kh_key(h, k), &ret);
                if (ret) { // absent
                    kh_key(h2, k2) = strdup(kh_key(h, k));
                    kh_val(h2, k2) = kh_val(h, k);
                } else kh_val(h2, k2) += kh_val(h, k);
            }
        }
    }
    k = kh_get(str, h, "foo_1"); printf("%d", kh_val(h, k));
    k = kh_get(str, h, "foo_9999"); printf(" %d", kh_val(h, k));
    k = kh_get(str, h2, "foo_1"); printf(" %d", kh_val(h2, k));
    k = kh_get(str, h2, "foo_9999"); printf(" %d\n", kh_val(h2, k));
    for (k = kh_begin(h); k != kh_end(h); ++k)
        if (kh_exist(h, k)) free((char*)kh_key(h, k));
    for (k = kh_begin(h2); k != kh_end(h2); ++k)
        if (kh_exist(h2, k)) free((char*)kh_key(h2, k));
    kh_destroy(str, h);
    kh_destroy(str, h2);
    return 0;
}*/