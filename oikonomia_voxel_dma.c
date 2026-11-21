
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdint.h>
#include <time.h>
#include <stdlib.h>

static inline uint64_t now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

int stream_voxel_batch_bar0(
    const uint8_t* blocks,
    size_t num_blocks,
    double* latencies_ns,
    uint32_t* hw_signatures
) {
    const char* pci_path = "/sys/bus/pci/devices/0000:34:00.0/resource0";
    int fd = open(pci_path, O_RDWR | O_SYNC);
    if (fd < 0) return -1;

    void* map_base = mmap(0, 4096, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (map_base == MAP_FAILED) {
        close(fd);
        return -2;
    }

    volatile uint32_t* reg = (volatile uint32_t*)map_base;
    for (int w = 0; w < 1000; w++) {
        reg[0] = 0xA5A5A5A5U;
    }

    for (size_t b = 0; b < num_blocks; b++) {
        const uint32_t* word_ptr = (const uint32_t*)(blocks + b * 4096ULL);
        uint64_t t0 = now_ns();
        for (int i = 0; i < 1024; i++) {
            reg[0] = word_ptr[i];
        }
        reg[4] = 0x00000001U;
        uint32_t out_val = reg[8];
        uint64_t t1 = now_ns();

        latencies_ns[b] = (double)(t1 - t0);
        hw_signatures[b] = out_val;
    }

    munmap(map_base, 4096);
    close(fd);
    return 0;
}
