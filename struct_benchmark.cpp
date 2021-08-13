#include "stdio.h"
#include "omp.h"
//#define USE_GPU
#include <ftrace.h>
#include <chrono>
#include <cmath>
#define ITER 10

#define USE_NEC

struct vec {
    float x, y, z;
};


template<typename T>
void test_combined_bw_nec(T *data1, T *data2, size_t length) {
    std::chrono::high_resolution_clock::time_point tstart;
    std::chrono::high_resolution_clock::time_point tstop;
    std::chrono::duration<double, std::milli> time_span;
    double dt;


    //    tstart = std::chrono::high_resolution_clock::now();
    //(void) ftrace_region_begin("benchmark");
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    //op1 6 * sizeof(float) access
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}
//    tstop = std::chrono::high_resolution_clock::now();
//    time_span = tstop - tstart;
//    dt = time_span.count();
//    printf("op1 BW is %lf GB/s \n", 9 * 1000 * length * sizeof(float) / dt);
//
//    tstart = std::chrono::high_resolution_clock::now();
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    //op2 6 * sizeof(float) accessS
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}
//    tstop = std::chrono::high_resolution_clock::now();
//    time_span = tstop - tstart;
//    dt = time_span.count();
//    //printf("op2 BW is %lf GB/s \n", 6 * 1000 * length * sizeof(float) / dt);
//
//    tstart = std::chrono::high_resolution_clock::now();
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {

    //op3 6 * sizeof(float) access
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}
//    tstop = std::chrono::high_resolution_clock::now();
//    time_span = tstop - tstart;
//    dt = time_span.count();
//printf("op3 BW is %lf GB/s \n", 6 * 1000 * length * sizeof(float) / dt);

//    tstart = std::chrono::high_resolution_clock::now();
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    //op4 9 * sizeof(float) access
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}
//    tstop = std::chrono::high_resolution_clock::now();
//    time_span = tstop - tstart;
//    dt = time_span.count();
// printf("op4 BW is %lf GB/s \n", 6 * 1000 * length * sizeof(float) / dt);

//tstart = std::chrono::high_resolution_clock::now();
#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
for (int i = 0; i < length; i++) {
    //op5 9 * sizeof(float) access
    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
}
/*tstop = std::chrono::high_resolution_clock::now();
time_span = tstop - tstart;
dt = time_span.count();
//printf("op5 BW is %lf GB/s \n", 6 * 1000 * length * sizeof(float) / dt);

 */
//    (void) ftrace_region_end("benchmark");
}


//template<typename T>
//void test_x_bw_nec(T *data1, T *data2, size_t length) {
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
//    for (int i = 0; i < length; i++) {
//        data1[i].x = data2[i].x;
//    }
//}
//
//template<typename T>
//void test_xy_bw_nec(T *data1, T *data2, size_t length) {
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
//    for (int i = 0; i < length; i++) {
//        data1[i].x = data2[i].x;
//        data1[i].y = data2[i].y;
//    }
//}
//
//template<typename T>
//void test_xyz_bw_nec(T *data1, T *data2, size_t length) {
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
//    for (int i = 0; i < length; i++) {
//        data1[i].x = data2[i].x;
//        data1[i].y = data2[i].y;
//        data1[i].z = data2[i].z;
//    }
//}
//
//template<typename T>
//void test_bw_nec(T *data1, T *data2, size_t length) {
//#pragma _NEC ivdep
//#pragma _NEC cncall
//#pragma omp parallel for
//    for (int i = 0; i < length; i++) {
//        data1[i] = data2[i];
//    }
//}

template<typename T>
void allocData(T **data, size_t length) {
#ifdef USE_NEC
    *data = (T *) aligned_alloc(sizeof(T), length * sizeof(T));
#endif
#ifdef USE_GPU
    SAFE_CALL(cudaMalloc((void**)data, length * sizeof(T)));
#endif
}

template<typename T>
void freeData(T *data) {
#ifdef USE_NEC
    free(data);
#endif
#ifdef USE_GPU
    cudaFree((void*)data);
#endif
}

int main(int argc, char** argv) {
    int coef = atoi(argv[1]);
    size_t length = 32 * 32 * 32 * coef;
    //    printf("Current length is %ld  ", length);

    float time_copy;

    //    float* dev_data1;
    //    float* dev_data2;
    //
    //    allocData<float>(&dev_data1, length);
    //    allocData<float>(&dev_data2, length);
    vec *dev_data1;
    vec *dev_data2;

    allocData<vec>(&dev_data1, length);
    allocData<vec>(&dev_data2, length);

#pragma _NEC ivdep
#pragma _NEC cncall
#pragma omp parallel for
    for (int i = 0; i < length; i++) {
        dev_data1[i].x = 1.0f;
        dev_data2[i].x = 2.0f;
        dev_data1[i].y = 1.0f;
        dev_data2[i].y = 2.0f;
        dev_data1[i].z = 1.0f;
        dev_data2[i].z = 2.0f;
    }

#ifdef USE_GPU
SAFE_CALL(cudaMemcpy(dev_data1, data1, sizeof(vec) * length, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(dev_data2, data2, sizeof(vec) * length, cudaMemcpyHostToDevice));
    //    SAFE_CALL(cudaMemcpy(dev_data1, data1, sizeof(float) * length, cudaMemcpyHostToDevice));
    //    SAFE_CALL(cudaMemcpy(dev_data2, data2, sizeof(float) * length, cudaMemcpyHostToDevice));
#endif

int blockSize = 1024;
int numBlocks = (length + blockSize - 1) / blockSize;
double final_bw;

for (int j = 0; j < ITER; j++) {
#ifdef USE_NEC
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

    //time_copy = time_copy * 1000;

#endif

//printf("lenght is %ld\n", length);
//printf("time  is %lf\n", dt);
double bw = 5 * 3 * 1000 * length * sizeof(vec) / (dt * std::pow(10, 9));
if (j == ITER - 1) {
    final_bw = bw;
}
//printf(".x BW is %lf GB/s \n", bw);
}

printf("  %ld\t%lf  ", length, final_bw);

/*
begin = omp_get_wtime();
#ifdef USE_GPU
SAFE_KERNEL_CALL((test_xy_bw<<<numBlocks, blockSize>>>(dev_data1, dev_data2)));
#endif
#ifdef USE_NEC
test_xy_bw_nec(data1, data2, length);
#endif
end = omp_get_wtime();
printf(".x .y BW is %f GB/s \n ", 4 * length * sizeof(float)/(end - begin));


begin = omp_get_wtime();
#ifdef USE_GPU
SAFE_KERNEL_CALL((test_xyz_bw<<<numBlocks, blockSize>>>(dev_data1, dev_data2)));
#endif
#ifdef USE_NEC
test_xyz_bw_nec(data1, data2, length);
#endif
end = omp_get_wtime();
printf(".x .y .z BW is %f GB/s\n", 6 * length * sizeof(float)/(end - begin));
*/


#ifdef USE_GPU
cudaMemcpy(data1, dev_data1, sizeof(float) * length, cudaMemcpyDeviceToHost);
cudaMemcpy(data2, dev_data2, sizeof(float) * length, cudaMemcpyDeviceToHost);
#endif


printf("%f", dev_data1[102400].x);
printf("%f", dev_data1[34].y);
printf("%f\n", dev_data1[54].z);

freeData(dev_data1);
freeData(dev_data2);
//    delete[] data1;
//    delete[] data2;
}
