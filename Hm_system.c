///This project simulates the ability of Hamming code correction.
///The input file contains: the number of rows of the H matrix,
///the number of columns of the H matrix, the H matrix itself,
///the number of error patterns, and the error pattern itself.

///This project needs to convert the parity check matrix into a generator matrix first,
///then a set of sequences is generated representing the desired information.
///Multiply the information sequence by the G matrix to complete the encoding,
///The information sequence is multiplied by the G matrix to complete the encoding,
///and the error pattern is added to observe the error correction ability of the Hamming code.
///Multiply the codeword plus the error pattern by the H matrix to get its syndrome.
///Make error corrections according to syndrome.
///Using the difference in error patterns to observe the error correction ability of Hamming codes

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>


///This function is used to convert H matrix to G matrix
void HtoG(int *H, int H_n, int H_k, int* G, int G_n, int G_k) {

    ///Gaussian Elimination
    ///--------------------------------------------------------------

    //1. Make sure the first element of the array G_matrix is 1.
    if (H[(H_k - 1) * H_n + (H_n - 1)] == 0) {
        for (int i = H_k - 1; i >= 0 ; i++) {
            if (H[i * H_n + H_k] == 1) {
                for (int j = 0; j < H_n; j++) {
                    H[H_k * H_n + j] = H[H_k * H_n + j] ^ H[i * H_n + j];
                }
            break;
            }
        }
    }

    //2. Find the first 1 that appears in each row,
    //and then eliminate the 1 of the same column and different row
    for (int i = H_k - 1; i >= 0; i--) {
        for(int j = H_n - 1; j >= 0; j--) {

            if (H[i * H_n + j] == 1) {

                for (int a = H_k - 1; a >= 0; a--) {

                    if (a == i || H[a * H_n + j] == 0)
                        continue;

                    else if (H[a * H_n + j] == 1) {
                        for (int jj = H_n - 1; jj >= 0; jj--) {
                            H[a * H_n + jj] = H[a * H_n + jj] ^ H[i * H_n + jj];
                        }
                    }
                }
                break;
            }
        }
    }

    //3. Make sure G is stepped
    int* I_step = malloc(H_k * sizeof(int));
    int index = H_k - 1;
    for (int i = H_k - 1; i >= 0; i--) {
        for (int j = H_n - 1; j >= 0; j--) {
            if (H[i * H_n + j] == 1) {
                I_step[index] = j;
                index--;
                break;
            }
        }
    }

    int* temp_mat = malloc(H_n * sizeof(int));
    int temp;
    for (int i = 0; i < (H_k - 1); i++) {
        for (int j = 0; j < (H_k - 1 - i); j++) {
            if (I_step[j] > I_step[j+1]) {
                temp = I_step[j];
                I_step[j] = I_step[j+1];
                I_step[j+1] = temp;

                for (int x = 0; x < H_n; x++) {
                    temp_mat[x] = H[j * H_n + x];
                    H[j * H_n + x] = H[(j + 1) * H_n + x];
                    H[(j + 1) * H_n + x] = temp_mat[x];
                }
            }
        }
    }
    ///Complete Gaussian Elimination



    ///----------------------------------------------------------------------
    ///Transform H to G
    //Find the elements in G
    int* I_loc = malloc(H_n * sizeof(int));
    for (int j = 0; j < H_n; j++) {
            I_loc[j] = 0;
    }
    for (int i = H_k - 1; i >= 0; i--) {
        for (int j = H_n - 1; j >= 0; j--) {
            if (H[i * H_n + j] == 1) {
                I_loc[j] = 1;
                break;
            }
        }
    }

    //1. Complete the part of the identity matrix first
    int* I_mat = malloc(G_k * G_k * sizeof(int));
    for (int i = 0; i < G_k; i++) {
        for (int j = 0; j < G_k; j++) {
            if (i == j)
                I_mat[i * G_k + j] = 1;

            else
                I_mat[i * G_k + j] = 0;
        }
    }

    //2. The part of the non-identity matrix is obtained by transposing the corresponding elements of the H matrix.
    int* A_Tran = malloc(G_k * (G_n - G_k) * sizeof(int));
    int row_count = 0;
    int col_count = 0;
    for (int i = 0; i < H_k; i++) {

        row_count = 0;
        for (int j = 0; j < H_n; j++) {

            if (I_loc[j] == 0) {
                A_Tran[row_count * H_k + col_count] = H[i * H_n + j];
                row_count++;
            }
        }
        col_count++;
    }

    //Combine the matrices obtained in step 1 and step 2 to obtain the H matrix to complete the conversion.
    for (int i = 0; i < G_k; i++) {
        for (int j = 0; j < G_n; j++) {
            G[i * G_n + j] = 0;
        }
    }

    int count = 0;
    for (int i = 0; i < G_k; i++) {
        for (int  j= 0; j < G_n; j++) {
            if (I_loc[j] == 1) {
                G[i * G_n + j] = A_Tran[count];
                count++;
            }
        }
    }

    count = 0;
    for (int i = 0; i < G_k; i++) {
        for (int  j= 0; j < G_n; j++) {
            if (I_loc[j] == 0) {
                G[i * G_n + j] = I_mat[count];
                count++;
            }
        }
    }

    free(I_step);
    free(temp_mat);
    free(I_loc);
    free(I_mat);
    free(A_Tran);
}

///This function is used to generate information sequence
//m: How many error patterns to simulate
void u_gen (int m, int G_k, int* u) {
    u[0] = 1;
    u[1] = 0;
    u[2] = 0;
    u[3] = 0;
    u[4] = 0;
    u[5] = 0;

    for(int i = 6; i < m * G_k; i++) {
        u[i] = u[i - 6] ^ u[i - 5];
    }
}

///This function is used for matrix multiplication (Binary mode)
void matrix_multiple(int* m_1, int* m_2, int* result, int m_1_row, int m_2_row, int m_2_col) {

    for (int i = 0; i < m_1_row; i++) {
        for (int j = 0; j < m_2_col; j++) {

            result[i * m_2_col + j] = 0;
            for (int k = 0; k < m_2_row; k++) {
                result[i * m_2_col + j] = result[i * m_2_col + j] ^ (m_1[i * m_2_row + k] * m_2[k * m_2_col + j]);
            }
        }
    }

}

///This function is used for Hamming detector
//Check syndrome patterns to choose which bit needs to be flip.
void Hm_detector(int* y, int* syndrome, int* H, int H_k, int H_n, int m) {

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < H_n; j++) {
            for (int k = 0; k < H_k; k++) {
                if (syndrome[k * m + i] != H[k * H_n + j])
                    break;
                else if (syndrome[k * m + i] == H[k * H_n + j] && k != H_k - 1)
                    continue;
                else if (syndrome[k * m + i] == H[k * H_n + j] && k == H_k - 1)
                    y[i * H_n + j] = !y[i * H_n + j];
            }
        }
    }

}

const char* filename = "Sim.txt";

int main(void) {

    ///Read file & transform the contents from "char" to "int"
    ///--------------------------------------------------------------
    FILE *in_file = fopen(filename, "r");
    if (!in_file) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    struct stat sb;
    if (stat(filename, &sb) == -1) {
        perror("stat");
        exit(EXIT_FAILURE);
    }

    char *file_contents = malloc(sb.st_size);
    int int_data[500];
    int index = 0;

    // Scan all "char" in the file but ignore blank & the "chars" between "%" and newline character
    const char* d = " \n";
    char *p;
    while (fscanf(in_file, "%[^%]%*[^\n] ", file_contents) != EOF) {
        p = strtok(file_contents, d);
        //printf("%s\n", p);
        while (p != NULL) {
            //printf("%s ", p);
            int_data[index] = atoi(p);
            //printf("%d\n", int_data[index]);
            index++;
            p = strtok(NULL, d);
        }
    }
/*
    for(int i = 0; i < 52; i++) {
        printf("%d", int_data[i]);
    }
*/
    fclose(in_file);


    ///------------------------------------------------------------------------
    //Form the H matrix to store the int input transformed from file
    int H_n = int_data[0];
    int H_k = int_data[1];
    int *H = malloc(H_n * H_k * sizeof(int));

    index = 2;
    for (int i = 0; i < H_k; i++) {
        for (int j = 0; j < H_n; j++) {
            H[i * H_n + j] = int_data[index];
            index++;
        }
    }

    //m: How many error patterns to simulate
    int m = int_data[index];
    index++;

    //Dimension of G matrix
    int G_n = H_n;
    int G_k = H_n - H_k;

    //error patterns
    int *error = malloc(m * G_n * sizeof(int));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < G_n; j++) {
            error[i * G_n + j] = int_data[index];
            index++;
        }
    }

    ///Convert G matrix from H matrix
    int *G = malloc(G_k * G_n * sizeof(int));
    HtoG(H, H_n, H_k, G, G_n, G_k);

    ///Generate information sequence
    int *u = malloc(m * G_k * sizeof(int));
    u_gen(m, G_k, u);

    ///Generate codewords
    int *x = calloc(m * G_n, sizeof(int));
    matrix_multiple(u, G, x, m, G_k, G_n);

    ///Add error patterns to codewords
    int *y = malloc(m * G_n * sizeof(int));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < G_n; j++) {
            y[i * G_n + j] = x[i * G_n + j] ^ error[i * G_n + j];
        }
    }

    ///Transpose codewords plus error pattern to obtain syndromes
    int *y_trans = malloc(H_n * m * sizeof(int));
    for (int i = 0; i < H_n; i++) {
        for (int j = 0; j < m; j++) {
            y_trans[i * m + j] = y[j * H_n + i];
        }
    }
    int *syndrome = malloc(H_k * m * sizeof(int));
    matrix_multiple(H, y_trans, syndrome, H_k, H_n, m);

    ///Detect
    Hm_detector(y, syndrome, H, H_k, H_n, m);

    ///Get the decoded information sequences
    int *u_hat = malloc(m * G_k * sizeof(int));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < G_k; j++) {
            u_hat[i * G_k + j] = y[i * G_n + j];
        }
    }


    ///----------------------------------------------------------------------------
    ///Write the decoded information sequences into an output txt file
    char *num = malloc(m * sizeof(char));
    for (int i = 0; i < m; i++) {
        sprintf(&num[i], "%d", i);
    }

    char u_hat_char[500] = {};
    index = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < G_k; j++) {
            sprintf(&u_hat_char[index], "%d", u_hat[i * G_k + j]);
            index++;
            u_hat_char[index] = ' ';
            index++;
        }

        //Notes
        u_hat_char[index] = '%';
        index++;
        u_hat_char[index] = '%';
        index++;
        u_hat_char[index] = 'E';
        index++;
        u_hat_char[index] = 's';
        index++;
        u_hat_char[index] = 't';
        index++;
        u_hat_char[index] = '.';
        index++;
        u_hat_char[index] = ' ';
        index++;
        u_hat_char[index] = 'u';
        index++;
        u_hat_char[index] = num[i];
        index++;
        u_hat_char[index] = '\n';
        index++;
    }
/*
    for (int r = 0; r < 500; r++) {
        printf("%c", u_hat_char[r]);
    }
*/

    printf("Complete! Please see the u.txt in the folder.\n");

    FILE *out_txt;
    out_txt = fopen("u.txt","w");
    fprintf(out_txt, u_hat_char);
    fclose(out_txt);

    getchar();

    free(G);
    free(H);
    free(u);
    free(x);
    free(error);
    free(y);
    free(y_trans);
    free(syndrome);
    free(u_hat);
    free(num);


    free(file_contents);

    exit(EXIT_SUCCESS);
}
