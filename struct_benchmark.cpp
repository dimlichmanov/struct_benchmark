#include <cstdio>
#include "omp.h"

#include <ftrace.h>
#include <chrono>
#include <cmath>

#define ITER 10

//#define USE_FLOAT
#define USE_VEC

struct vec {
    float x, y, z;
};


template<typename T>
void test_combined_bw_nec(T *data1, T *data2, size_t length) {
    //(void) ftrace_region_begin("benchmark");
#ifdef USE_VEC
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }

//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }

//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        data1[i].x = (data1[i].x - data2[i].x) / 2;
        data1[i].y = (data1[i].y - data2[i].y) / 2;
        data1[i].z = (data1[i].z - data2[i].z) / 2;
    }
#endif

#ifdef USE_FLOAT
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }

//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
    }

//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
for (int i = 0; i < length; i++) {
        data1[i] = (data1[i]- data2[i]) / 2;
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
    size_t length = 32 * 8 * 1 * coef;

#ifdef USE_FLOAT
    float *dev_data1;
    float *dev_data2;

    allocData<float>(&dev_data1, length);
    allocData<float>(&dev_data2, length);
#endif
#ifdef USE_VEC
    vec *dev_data1;
    vec *dev_data2;

    allocData<vec>(&dev_data1, length);
    allocData<vec>(&dev_data2, length);
#endif

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        #ifdef USE_FLOAT
        dev_data1[i] = 1.0f;
        dev_data2[i] = 2.0f;
#endif
#ifdef USE_VEC
        dev_data1[i].x = 1.0f;
        dev_data2[i].x = 2.0f;
        dev_data1[i].y = 1.0f;
        dev_data2[i].y = 2.0f;
        dev_data1[i].z = 1.0f;
        dev_data2[i].z = 2.0f;
#endif
    }

    double final_bw;

    for (int j = 0; j < ITER; j++) {
        std::chrono::high_resolution_clock::time_point tstart;
        tstart = std::chrono::high_resolution_clock::now();
        //test_bw_nec(dev_data1, dev_data2, length);
        //test_x_bw_nec(dev_data1, dev_data2, length);
        //test_xy_bw_nec(dev_data1, dev_data2, length);
        test_combined_bw_nec(dev_data1, dev_data2, length);
        //test_xyz_bw_nec(dev_data1, dev_data2, length);
        std::chrono::high_resolution_clock::time_point tstop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_span = tstop - tstart;
        double dt = time_span.count(); // in milliseconds

#ifdef USE_FLOAT
        double bw = 5.0 * 3 * 1000 * length * sizeof(float) / (dt * std::pow(10, 9));
#endif
#ifdef USE_VEC
        double bw = 5.0 * 3 * 1000 * length * sizeof(vec) / (dt * std::pow(10, 9));
#endif
        if (j == ITER - 1) {
            final_bw = bw;
        }
//printf(".x BW is %lf GB/s \n", bw);
    }

    printf("  %ld\t%lf  ", length, final_bw);

#ifdef USE_FLOAT
    printf("%f", dev_data1[length - 1]);
    printf("%f", dev_data1[34]);
    printf("%f\n", dev_data1[54]);
#endif
#ifdef USE_VEC
    printf("%f", dev_data1[length - 1].x);
    printf("%f", dev_data1[34].y);
    printf("%f\n", dev_data1[54].z);
#endif

    freeData(dev_data1);
    freeData(dev_data2);
}
