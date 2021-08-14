#include <cstdio>
#include "omp.h"

#include <chrono>
#include <cmath>

#define ITER 10

//#define USE_NEC
#define USE_GPU
#define USE_FLOAT
//#define USE_VEC

#define SAFE_CALL( CallInstruction ) { \
cudaError_t cuerr = CallInstruction; \
cudaDeviceSynchronize();\
if(cuerr != cudaSuccess) { \
printf("CUDA error: %s at call \"" #CallInstruction "\"\n", cudaGetErrorString(cuerr)); \
throw "error in CUDA API function, aborting..."; \
} \
}

#define SAFE_KERNEL_CALL( KernelCallInstruction ){ \
KernelCallInstruction; \
cudaError_t cuerr = cudaGetLastError(); \
cudaDeviceSynchronize();\
if(cuerr != cudaSuccess) { \
printf("CUDA error in kernel launch: %s at kernel \"" #KernelCallInstruction "\"\n", cudaGetErrorString(cuerr)); \
throw "error in CUDA kernel launch, aborting..."; \
} \
if(cuerr != cudaSuccess) { \
printf("CUDA error in kernel execution: %s at kernel \"" #KernelCallInstruction "\"\n", cudaGetErrorString(cuerr)); \
throw "error in CUDA kernel execution, aborting..."; \
} \
}

struct vec {
    float x, y, z;
};


template<typename T>
__global__ void test_combined_bw_gpu(T *data1, T *data2, size_t length) {
    //(void) ftrace_region_begin("benchmark");
    int i = blockDim.x * blockIdx.x + threadIdx.x;

#ifdef USE_VEC

    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;

    //cudaDeviceSynchronize();

    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;

    //cudaDeviceSynchronize();

    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;

    //cudaDeviceSynchronize();

    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;

    //cudaDeviceSynchronize();

    data1[i].x = (data1[i].x - data2[i].x) / 2;
    data1[i].y = (data1[i].y - data2[i].y) / 2;
    data1[i].z = (data1[i].z - data2[i].z) / 2;
#endif

#ifdef USE_FLOAT

    data1[i] = (data1[i]- data2[i]) / 2;
    //cudaDeviceSynchronize();

    data1[i] = (data1[i]- data2[i]) / 2;
    //cudaDeviceSynchronize();

    data1[i] = (data1[i]- data2[i]) / 2;
    //cudaDeviceSynchronize();

    data1[i] = (data1[i]- data2[i]) / 2;
    //cudaDeviceSynchronize();

    data1[i] = (data1[i]- data2[i]) / 2;
    //cudaDeviceSynchronize();
#endif
//    (void) ftrace_region_end("benchmark");
}


template<typename T>
__global__ void test_x_bw(T *data1, T *data2) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    data1[i].x = data2[i].x;
}

template<typename T>
__global__ void test_xy_bw(T *data1, T *data2) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    data1[i].x = data2[i].x;
    data1[i].y = data2[i].y;
}

template<typename T>
__global__ void test_xyz_bw(T *data1, T *data2) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    data1[i].x = data2[i].x;
    data1[i].y = data2[i].y;
    data1[i].z = data2[i].z;
}

template<typename T>
void allocData(T **data, size_t length) {
#ifdef USE_NEC
    *data = (T *) aligned_alloc(sizeof(T), length * sizeof(T));
#endif
#ifdef USE_GPU
    SAFE_CALL(cudaMalloc((void **) data, length * sizeof(T)));
#endif
}

template<typename T>
void freeData(T *data) {
#ifdef USE_NEC
    free(data);
#endif
#ifdef USE_GPU
    cudaFree((void *) data);
#endif
}

int main(int argc, char **argv) {
    int coef = atoi(argv[1]);
    size_t length = 32 * 32 * 32 * coef;

#ifdef USE_FLOAT
    float *data1 = new float[length];
    float *data2 = new float[length];
    float *dev_data1;
    float *dev_data2;

    allocData<float>(&dev_data1, length);
    allocData<float>(&dev_data2, length);
#endif
#ifdef USE_VEC
    vec *data1 = new vec[length];
    vec *data2 = new vec[length];
    vec *dev_data1;
    vec *dev_data2;

    allocData<vec>(&dev_data1, length);
    allocData<vec>(&dev_data2, length);
#endif

#pragma omp parallel for
    for (int i = 0; i < length; i++) {
#ifdef USE_FLOAT
        data1[i] = 1.0f;
        data2[i] = 2.0f;
#endif
#ifdef USE_VEC
        data1[i].x = 1.0f;
        data2[i].x = 2.0f;
        data1[i].y = 1.0f;
        data2[i].y = 2.0f;
        data1[i].z = 1.0f;
        data2[i].z = 2.0f;
#endif
    }

#ifdef USE_FLOAT
    SAFE_CALL(cudaMemcpy(dev_data1, data1, sizeof(float) * length, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(dev_data2, data2, sizeof(float) * length, cudaMemcpyHostToDevice));
#endif
#ifdef USE_VEC
    SAFE_CALL(cudaMemcpy(dev_data1, data1, sizeof(vec) * length, cudaMemcpyHostToDevice));
    SAFE_CALL(cudaMemcpy(dev_data2, data2, sizeof(vec) * length, cudaMemcpyHostToDevice));
#endif

    int blockSize = 1024;
    int numBlocks = (length + blockSize - 1) / blockSize;
    double final_bw;

    for (int j = 0; j < ITER; j++) {
        std::chrono::high_resolution_clock::time_point tstart;
        tstart = std::chrono::high_resolution_clock::now();
        //test_bw_nec(dev_data1, dev_data2, length);
        //test_x_bw_nec(dev_data1, dev_data2, length);
        //test_xy_bw_nec(dev_data1, dev_data2, length);
        SAFE_KERNEL_CALL((test_combined_bw_gpu<<<numBlocks, blockSize>>>(dev_data1, dev_data2, length)));
        //test_xyz_bw_nec(dev_data1, dev_data2, length);
        std::chrono::high_resolution_clock::time_point tstop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_span = tstop - tstart;
        double dt = time_span.count(); // in milliseconds

#ifdef USE_FLOAT
    double bw = 5 * 3 * 1000 * length * sizeof(float) / (dt * std::pow(10, 9));
#endif
#ifdef USE_VEC
    double bw = 5 * 3 * 1000 * length * sizeof(vec) / (dt * std::pow(10, 9));
#endif
        if (j == ITER - 1) {
            final_bw = bw;
        }
    }

    printf("  %ld\t%lf  ", length, final_bw);


#ifdef USE_FLOAT
    SAFE_CALL(cudaMemcpy(data1, dev_data1, sizeof(float) * length, cudaMemcpyDeviceToHost));
    SAFE_CALL(cudaMemcpy(data2, dev_data2, sizeof(float) * length, cudaMemcpyDeviceToHost));
#endif
#ifdef USE_VEC
    SAFE_CALL(cudaMemcpy(data1, dev_data1, sizeof(vec) * length, cudaMemcpyDeviceToHost));
    SAFE_CALL(cudaMemcpy(data2, dev_data2, sizeof(vec) * length, cudaMemcpyDeviceToHost));
#endif

#ifdef USE_FLOAT
    printf("%f", data1[102400]);
    printf("%f", data1[34]);
    printf("%f\n", data1[54]);
#endif
    #ifdef USE_VEC
    printf("%f", data1[102400].x);
    printf("%f", data1[34].y);
    printf("%f\n", data1[54].z);
    #endif

    freeData(dev_data1);
    freeData(dev_data2);
    delete[] data1;
    delete[] data2;
}
