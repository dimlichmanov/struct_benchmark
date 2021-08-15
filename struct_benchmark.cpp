#include <cstdio>
#include "omp.h"

#include <ftrace.h>
#include <chrono>
#include <cmath>

#define ITER 10

#define USE_FLOAT
//#define USE_VEC

struct vec {
    float x, y, z;
};

template<typename T>
void test_combined_bw_nec(T *data1, T *data2, size_t length) {
    //(void) ftrace_region_begin("benchmark");
#ifdef USE_VEC
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }
#endif

#ifdef USE_FLOAT
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }
#endif
//    (void) ftrace_region_end("benchmark");
}

template<typename T>
void test_complex_bw_nec(T *data1, T *data2, T *data3, T *data4, T *data5, T *data6, T *data7, T *data8, T *data9, T *data10, size_t length) {
    //(void) ftrace_region_begin("benchmark");
#ifdef USE_VEC
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}
#endif

#ifdef USE_FLOAT
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {

#pragma _NEC unroll_completely
for (size_t nx = 0; nx < 256; nx++) {
        data1[i] = (data1[i] + data2[i]) / 2;
        data2[i] = (data2[i] + data3[i]) / 2;
        data3[i] = (data3[i] + data4[i]) / 2;
        data4[i] = (data4[i] + data5[i]) / 2;
        data5[i] = (data5[i] + data6[i]) / 2;
        data6[i] = (data6[i] + data7[i]) / 2;
        data7[i] = (data7[i] + data8[i]) / 2;
        data8[i] = (data8[i]- data9[i]) / 2;
        data9[i] = (data9[i]- data10[i]) / 2;
        data10[i] = (data1[i] + data10[i]) / 2;
    }
}

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
#pragma _NEC unroll_completely
    for (size_t nx = 0; nx < 256; nx++) {
        data1[i] = (data1[i] + data2[i]) / 2;
        data2[i] = (data2[i] + data3[i]) / 2;
        data3[i] = (data3[i] + data4[i]) / 2;
        data4[i] = (data4[i] + data5[i]) / 2;
        data5[i] = (data5[i] + data6[i]) / 2;
        data6[i] = (data6[i] + data7[i]) / 2;
        data7[i] = (data7[i] + data8[i]) / 2;
        data8[i] = (data8[i]- data9[i]) / 2;
        data9[i] = (data9[i]- data10[i]) / 2;
        data10[i] = (data1[i] + data10[i]) / 2;
    }
}
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
#pragma _NEC unroll_completely
    for (size_t nx = 0; nx < 256; nx++) {
        data1[i] = (data1[i] + data2[i]) / 2;
        data2[i] = (data2[i] + data3[i]) / 2;
        data3[i] = (data3[i] + data4[i]) / 2;
        data4[i] = (data4[i] + data5[i]) / 2;
        data5[i] = (data5[i] + data6[i]) / 2;
        data6[i] = (data6[i] + data7[i]) / 2;
        data7[i] = (data7[i] + data8[i]) / 2;
        data8[i] = (data8[i]- data9[i]) / 2;
        data9[i] = (data9[i]- data10[i]) / 2;
        data10[i] = (data1[i] + data10[i]) / 2;
   }
}
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
#pragma _NEC unroll_completely
    for (size_t nx = 0; nx < 256; nx++) {
        data1[i] = (data1[i] + data2[i]) / 2;
        data2[i] = (data2[i] + data3[i]) / 2;
        data3[i] = (data3[i] + data4[i]) / 2;
        data4[i] = (data4[i] + data5[i]) / 2;
        data5[i] = (data5[i] + data6[i]) / 2;
        data6[i] = (data6[i] + data7[i]) / 2;
        data7[i] = (data7[i] + data8[i]) / 2;
        data8[i] = (data8[i]- data9[i]) / 2;
        data9[i] = (data9[i]- data10[i]) / 2;
        data10[i] = (data1[i] + data10[i]) / 2;
    }
}

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
  #pragma _NEC unroll_completely
    for (size_t nx = 0; nx < 256; nx++) {
        data1[i] = (data1[i] + data2[i]) / 2;
        data2[i] = (data2[i] + data3[i]) / 2;
        data3[i] = (data3[i] + data4[i]) / 2;
        data4[i] = (data4[i] + data5[i]) / 2;
        data5[i] = (data5[i] + data6[i]) / 2;
        data6[i] = (data6[i] + data7[i]) / 2;
        data7[i] = (data7[i] + data8[i]) / 2;
        data8[i] = (data8[i]- data9[i]) / 2;
        data9[i] = (data9[i]- data10[i]) / 2;
        data10[i] = (data1[i] + data10[i]) / 2;
    }
}
#endif
//    (void) ftrace_region_end("benchmark");
}


template<typename T>
void test_x_bw_nec(T *data1, T *data2, size_t length) {
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = data2[i].x;
    }
}

template<typename T>
void test_xy_bw_nec(T *data1, T *data2, size_t length) {
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = data2[i].x;
        data1[i].y = data2[i].y;
    }
}

template<typename T>
void test_xyz_bw_nec(T *data1, T *data2, size_t length) {
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = data2[i].x;
        data1[i].y = data2[i].y;
        data1[i].z = data2[i].z;
    }
}

template<typename T>
void allocData(T **data, size_t length) {
    *data = (T *) aligned_alloc(sizeof(T), length * sizeof(T));
}

template<typename T>
void freeData(T *data) {
    free(data);
}

int main(int argc, char **argv) {
    int coef = atoi(argv[1]);
    size_t length = 32 * 32 * 32 * coef;

#ifdef USE_FLOAT
    float *dev_data1;
    float *dev_data2;
    float *dev_data3;
    float *dev_data4;
    float *dev_data5;
    float *dev_data6;
    float *dev_data7;
    float *dev_data8;
    float *dev_data9;
    float *dev_data10;

    allocData<float>(&dev_data1, length);
    allocData<float>(&dev_data2, length);
    allocData<float>(&dev_data3, length);
    allocData<float>(&dev_data4, length);
    allocData<float>(&dev_data5, length);
    allocData<float>(&dev_data6, length);
    allocData<float>(&dev_data7, length);
    allocData<float>(&dev_data8, length);
    allocData<float>(&dev_data9, length);
    allocData<float>(&dev_data10, length);

#endif
#ifdef USE_VEC
    vec *dev_data1;
    vec *dev_data2;
    vec *dev_data3;
    vec *dev_data4;
    vec *dev_data5;
    vec *dev_data6;
    vec *dev_data7;
    vec *dev_data8;
    vec *dev_data9;
    vec *dev_data10;

    allocData<vec>(&dev_data1, length);
    allocData<vec>(&dev_data2, length);
    allocData<vec>(&dev_data3, length);
    allocData<vec>(&dev_data4, length);
    allocData<vec>(&dev_data5, length);
    allocData<vec>(&dev_data6, length);
    allocData<vec>(&dev_data7, length);
    allocData<vec>(&dev_data8, length);
    allocData<vec>(&dev_data9, length);
    allocData<vec>(&dev_data10, length);
#endif

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        #ifdef USE_FLOAT
        dev_data1[i] = 1.0f;
        dev_data2[i] = 2.0f;
        dev_data3[i] = 2.0f;
        dev_data4[i] = 2.0f;
        dev_data5[i] = 2.0f;
        dev_data6[i] = 2.0f;
        dev_data7[i] = 2.0f;
        dev_data8[i] = 2.0f;
        dev_data9[i] = 2.0f;
        dev_data10[i] = 2.0f;

#endif
#ifdef USE_VEC
        dev_data1[i].x = 1.0f;
        dev_data2[i].x = 2.0f;
        dev_data1[i].y = 1.0f;
        dev_data2[i].y = 2.0f;
        dev_data1[i].z = 1.0f;
        dev_data2[i].z = 2.0f;
        dev_data3[i].x = 1.0f;
        dev_data4[i].x = 2.0f;
        dev_data3[i].y = 1.0f;
        dev_data4[i].y = 2.0f;
        dev_data3[i].z = 1.0f;
        dev_data4[i].z = 2.0f;
        dev_data5[i].x = 1.0f;
        dev_data6[i].x = 2.0f;
        dev_data5[i].y = 1.0f;
        dev_data6[i].y = 2.0f;
        dev_data5[i].z = 1.0f;
        dev_data6[i].z = 2.0f;
        dev_data7[i].x = 1.0f;
        dev_data8[i].x = 2.0f;
        dev_data7[i].y = 1.0f;
        dev_data8[i].y = 2.0f;
        dev_data7[i].z = 1.0f;
        dev_data8[i].z = 2.0f;
        dev_data9[i].x = 1.0f;
        dev_data10[i].x = 2.0f;
        dev_data9[i].y = 1.0f;
        dev_data10[i].y = 2.0f;
        dev_data9[i].z = 1.0f;
        dev_data10[i].z = 2.0f;
#endif
    }

    double final_bw;

    for (int j = 0; j < ITER; j++) {
        std::chrono::high_resolution_clock::time_point tstart;
        tstart = std::chrono::high_resolution_clock::now();
        //test_bw_nec(dev_data1, dev_data2, length);
        //test_x_bw_nec(dev_data1, dev_data2, length);
        //test_xy_bw_nec(dev_data1, dev_data2, length);
        //test_combined_bw_nec(dev_data1, dev_data2, length);
        test_complex_bw_nec(dev_data1, dev_data2, dev_data3, dev_data4, dev_data5, dev_data6, dev_data7, dev_data8, dev_data9, dev_data10, length);
        //test_xyz_bw_nec(dev_data1, dev_data2, length);
        std::chrono::high_resolution_clock::time_point tstop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_span = tstop - tstart;
        double dt = time_span.count(); // in milliseconds


        /* 5 - num of sequential iterations
         * 3 - every instruction
         * 10 - number of rows(arrays)
         * 1000 - for msec
         * 256 - inner cycle (linearized)
        */

#ifdef USE_FLOAT
        double bw = 5.0 * 3 * 10 * 1000 * 256 * length * sizeof(float) / (dt * std::pow(10, 9));
#endif
#ifdef USE_VEC
        double bw = 5.0 * 3 * 10 * 1000 * 256 * length * sizeof(vec) / (dt * std::pow(10, 9));
#endif
        if (j == ITER - 1) {
            final_bw = bw;
        }
        printf(".x BW is %lf GB/s \n", bw);
    }

    printf("  %ld\t%lf  ", length, final_bw);

#ifdef USE_FLOAT
    printf("%f", dev_data1[length - 1]);
    printf("%f", dev_data2[34]);
    printf("%f\n", dev_data3[54]);
    printf("%f\n", dev_data4[54]);
    printf("%f\n", dev_data5[54]);
    printf("%f\n", dev_data6[54]);
    printf("%f\n", dev_data7[54]);
    printf("%f\n", dev_data8[54]);
    printf("%f\n", dev_data9[54]);

    printf("%f\n", dev_data10[54]);
#endif
#ifdef USE_VEC
    printf("%f", dev_data1[length - 1].x);
    printf("%f", dev_data2[34].y);
    printf("%f\n", dev_data3[54].z);
    printf("%f\n", dev_data4[54].z);
    printf("%f\n", dev_data5[54].z);
    printf("%f\n", dev_data6[54].z);
    printf("%f\n", dev_data7[54].z);
    printf("%f\n", dev_data8[54].z);
    printf("%f\n", dev_data9[54].z);
    printf("%f\n", dev_data10[54].z);
#endif

    freeData(dev_data1);
    freeData(dev_data2);
    freeData(dev_data3);
    freeData(dev_data4);
    freeData(dev_data5);
    freeData(dev_data6);
    freeData(dev_data7);
    freeData(dev_data8);
    freeData(dev_data9);
    freeData(dev_data10);
}
