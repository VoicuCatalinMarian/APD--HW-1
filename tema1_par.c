// VOICU CATALIN MARIAN 332 CC Tema1 APD

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

// Structura ce contine datele necesare pentru thread-uri
typedef struct{
    ppm_image *image;
    ppm_image *new_image;
    ppm_image **contour_map;
    unsigned char **grid;
    uint8_t sample[3];
    int step_x;
    int step_y;
    int P; // Numarul de thread-uri
    int id; // Id-ul thread-ului
    int check; // Verifica daca e nevoie de rescale
    pthread_barrier_t *barrier; // Bariera
} ThreadData;

// Functie ce returneaza minimul dintre doua numere
int min(int a, int b){
    if(a < b){
        return a;
    }
    return b;
}

ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

// Functie pentru paraleliazare
void *thread_function(void* arg){
    ThreadData *thread = (ThreadData*)arg;


    // Rescale_Image 


    // Rescale daca este nevoie folosind Bicubic Interpolation
    if(thread->check == 1)
    {
        // Calcularea intervalului de prelucrare pentru fiecare thread
        int start_rescale = thread->id * (double)thread->new_image->x / thread->P;
        int end_rescale = min((thread->id + 1) * (double)thread->new_image->x / thread->P, thread->new_image->x);
        
        for (int i = start_rescale; i < end_rescale; i++) 
        {
            for (int j = 0; j < thread->new_image->y; j++) 
            {
                float u = (float)i / (float)(thread->new_image->x - 1);
                float v = (float)j / (float)(thread->new_image->y - 1);
                sample_bicubic(thread->image, u, v, thread->sample);

                thread->new_image->data[i * thread->new_image->y + j].red = thread->sample[0];
                thread->new_image->data[i * thread->new_image->y + j].green = thread->sample[1];
                thread->new_image->data[i * thread->new_image->y + j].blue = thread->sample[2];
            }
        }
    }

    // Se asteapta ca toate thread-urile sa termine executia
    pthread_barrier_wait(thread->barrier); 


    // Sample_Grid


    int p = thread->new_image->x / thread->step_x; // numarul de linii
    int q = thread->new_image->y / thread->step_y; // numarul de coloane

    // Calcularea intervalului de prelucrare pentru fiecare thread
    int start_sample = thread->id * (double)p / thread->P;
    int end_sample = min((thread->id + 1) * (double)p / thread->P, p);
 
    for (int i = start_sample; i < end_sample; i++) 
    {
        for (int j = 0; j < q; j++) 
        {
            ppm_pixel curr_pixel = thread->new_image->data[i * thread->step_x * thread->new_image->y + j * thread->step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                thread->grid[i][j] = 0;
            } else {
                thread->grid[i][j] = 1;
            }
        }
    }

    // Se asteapta ca toate thread-urile sa termine executia
    pthread_barrier_wait(thread->barrier);

    thread->grid[p][q] = 0;

    // Calcularea intervalului de prelucrare pentru fiecare thread
    start_sample = thread->id * (double)p / thread->P;
    end_sample = min((thread->id + 1) * (double)p / thread->P, p);

    for (int i = start_sample; i < end_sample; i++) 
    {
        ppm_pixel curr_pixel = thread->new_image->data[i * thread->step_x * thread->new_image->y + thread->new_image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            thread->grid[i][q] = 0;
        } else {
            thread->grid[i][q] = 1;
        }
    }

    // Se asteapta ca toate thread-urile sa termine executia
    pthread_barrier_wait(thread->barrier);

    // Calcularea intervalului de prelucrare pentru fiecare thread
    start_sample = thread->id * (double)q / thread->P;
    end_sample = min((thread->id + 1) * (double)q / thread->P, q);

    for (int j = start_sample; j < end_sample; j++) 
    {
        ppm_pixel curr_pixel = thread->new_image->data[(thread->new_image->x - 1) * thread->new_image->y + j * thread->step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            thread->grid[p][j] = 0;
        } else {
            thread->grid[p][j] = 1;
        }
    }
    
    // Se asteapta ca toate thread-urile sa termine executia
    pthread_barrier_wait(thread->barrier);
    
    
    // March


    // Calcularea intervalului de prelucrare pentru fiecare thread
    int start_march = thread->id * (double)p / thread->P;
    int end_march = min((thread->id + 1) * (double)p / thread->P, p);

    for (int i = start_march; i < end_march; i++) 
    {
        for (int j = 0; j < q; j++) 
        {
            unsigned char k = 8 * thread->grid[i][j] + 4 * thread->grid[i][j + 1] + 2 * thread->grid[i + 1][j + 1] + 1 * thread->grid[i + 1][j];
            update_image(thread->new_image, thread->contour_map[k], i * thread->step_x, j * thread->step_y);
        }
    }

    // Se asteapta ca toate thread-urile sa termine executia
    pthread_barrier_wait(thread->barrier);

    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    int check = 0; // Verifica daca e nevoie de rescale
    int P = atoi(argv[3]); // Numarul de thread-uri

    ppm_image *image = read_ppm(argv[1]);
    ppm_image *new_image;

    int step_x = STEP; 
    int step_y = STEP;

    pthread_t threads[P]; // Initializarea thread-urilor
    ThreadData thread_structure[P]; // Initializarea structurii de date pentru thread-uri
    
    // Initializarea barierei
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, P);

    // Initializarea hartii de contur
    ppm_image **contour_map = init_contour_map();

    // Verific daca e nevoie de rescale
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        new_image = image;
    }
    else {
        // Alocarea memoriei pentru noua imagine
        new_image = (ppm_image *)malloc(sizeof(ppm_image));
        if (!new_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }

        new_image->x = RESCALE_X;
        new_image->y = RESCALE_Y;

        // Alocarea memoriei pentru pixeli
        new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
        if (!new_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }

        check = 1; // Daca e nevoie de rescale, check devine 1
    }

    int p = new_image->x / step_x; // Numarul de linii
    int q = new_image->y / step_y; // Numarul de coloane

    // Alocarea memoriei pentru grid
    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    // Alocarea memoriei pentru fiecare linie din grid cu numarul de coloane
    for(int i = 0 ; i <= p ; i++){
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // Initializarea structurii de date pentru fiecare thread
    for(int id = 0 ; id < P ; id++){
        thread_structure[id].image = image;
        thread_structure[id].new_image = new_image;
        thread_structure[id].contour_map = contour_map;
        thread_structure[id].grid = grid;
        thread_structure[id].P = P;
        thread_structure[id].id = id;
        thread_structure[id].check = check;
        thread_structure[id].step_x = step_x;
        thread_structure[id].step_y = step_y;
        thread_structure[id].barrier = &barrier;
    }

    // Crearea thread-urilor
    for (int i = 0; i < P; i++) {
        pthread_create(&(threads[i]), NULL, thread_function, &(thread_structure[i]));
    }

    // Asteptarea terminarii thread-urilor
    for (int i = 0; i < P; i++) {
        pthread_join(threads[i], NULL);
    }

    // Scrierea in fisier a noii imagini
    write_ppm(new_image, argv[2]);

    // Distrugerea barierei
    pthread_barrier_destroy(&barrier);

    // Eliberarea memoriei
    free_resources(new_image, contour_map, grid, step_x);

    return 0;
}