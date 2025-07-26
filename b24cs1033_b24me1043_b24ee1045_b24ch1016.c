#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX 10

void addition(int*ptr1,int*ptr2,int*ptr3,int m,int n,int p,int q,FILE *file3) 
{
    if(m==p&&n==q) //condition for addition of matrices
    {
        for(int k=0;k<m;k++) 
        {
            for(int l=0;l<n;l++) 
            {
                *(ptr3+k*n+l) = *(ptr1+k*n+l) + *(ptr2+k*n+l);
            }
        }
        fprintf(file3, "Resulting Matrix:\n");
        for (int r=0;r<m;r++) 
        {
            for (int s=0;s<n;s++) 
            {
                fprintf(file3,"%d ",*(ptr3+r*n+s)); // printing the matrix
            }
            fprintf(file3,"\n");
        }
    }
}

void subtraction(int *ptr1,int *ptr2,int *ptr3,int m,int n,int p,int q,FILE *file3) 
{
    if(m==p&&n==q) // condition for subtraction of matrices
    {
        for(int k=0;k<m;k++) 
        {
            for(int l=0;l<n;l++) 
            {
                *(ptr3+k*n+l)=*(ptr1+k*n+l)-*(ptr2+k*n+l);
            }
        }
        fprintf(file3,"Resulting Matrix:\n");
        for(int r=0;r<m;r++) 
        {
            for(int s=0;s<n;s++) 
            {
                fprintf(file3,"%d ", *(ptr3+r*n+s)); // printing the matrix
            }
            fprintf(file3,"\n");
        }
    }    
}

void multiplication(int *ptr1,int *ptr2,int *ptr3,int m,int n,int p,int q,FILE *file3)
{
    if (n!=p) // condition for multiplication of matrices
    {
        printf("Multiplication is not possible\n");
    }
    else
    {
        for(int i=0;i<m;i++) 
        {
            for(int j=0;j<q;j++) 
            {
                *(ptr3+i*q+j)= 0;
            }
        }

        for(int i=0;i<m;i++) 
        {
            for (int j=0;j<q;j++) 
            {
                for (int k=0;k<n;k++) 
                {
                    (*(ptr3+i*q+j)) += (*(ptr1+i*n+k)) * (*(ptr2+k*q+j));
                }
            }
        }

        fprintf(file3,"Resulting Matrix:\n");
        for (int i=0;i<m;i++) 
        {
            for (int j=0;j<q;j++) 
            {
                fprintf(file3,"%d",*(ptr3 + i * q + j)); // printing the result
            }
            fprintf(file3,"\n");
        }
    }
}

int determinant(int *arr,int n)
{
    int det=0;
    int cofactor[MAX][MAX];

    if(n==1) 
    {
        return arr[0];
    }

    if(n==2) 
    {
        return arr[0]*arr[3]-arr[1]*arr[2];
    }

    for(int c=0;c<n;c++) 
    {
        int subi=0;
        for(int i=1;i<n;i++) 
        {
            int subj=0;
            for(int j=0;j<n;j++) 
            {
                if(j==c) 
                continue;
                cofactor[subi][subj] = arr[i*n+j];
                subj++;
            }
            subi++;
        }
        det += (c%2 == 0 ? 1 : -1) * arr[c] * determinant((int *)cofactor,n - 1);
    }
    return det;
}

void transposingthematrix(int *matrix,int *transpose1,int m,int n)
{
    for(int i=0;i<m;i++) 
    {
        for(int j=0;j<n;j++) 
        {
            *((int *)transpose1 +j*m +i) = *((int *)matrix +i*n+j); // transpose logic
        }
    }
}

void traceofmatrix(int *A,int n,FILE *file3)
{
    int trace=0;
    for(int i=0; i<n; i++) 
    {
        trace+= A[i*n+i]; // adding diagonal elements
    }
    fprintf(file3,"Trace of matrix: %d\n",trace);
}

int DET11(int *matrix,int n)
{
    int DET1=0;

    if(n==1) 
    {
        return matrix[0];
    }

    if(n==2) 
    {
        return matrix[0]*matrix[3]-matrix[1]*matrix[2];
    }

    for(int q=0;q<n;q++) 
    {
        int submatrix[MAX][MAX];
        int subr=0;
        for(int i=1;i<n;i++) 
        {
            int subc=0;
            for(int j=0;j<n;j++) 
            {
                if(j==q) 
                continue;
                submatrix[subr][subc]=matrix[i*n+j];
                subc++;
            }
            subr++;
        }
        DET1 +=(q%2 == 0 ? 1 : -1) * matrix[q] * DET11((int *)submatrix,n - 1);
    }
    return DET1;
}

void findcofactor(int *matrix,int *cofactor,int n)
{
    for(int k=0; k<n;k++) 
    {
        for(int l=0;l<n;l++) 
        {
            int submatrix[MAX][MAX];
            int subr=0;
            
            for(int g=0;g<n;g++) 
            {
                if(g==k) 
                continue;
                int subc=0;
                for(int h=0;h<n;h++) 
                {
                    if(h==l) 
                    continue;
                    submatrix[subr][subc]=matrix[g*n+h];
                    subc++;
                }
                subr++;
            }
            cofactor[k*n+l]=((k+l)%2 == 0 ? 1 : -1) * DET11((int *)submatrix,n - 1);
        }
    }
}

void findadjoint(int *cofactor,int *adjoint,int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            *(adjoint+i*n+j)= *(cofactor+j*n+i);
        }
    }
}

int findinverse(int *matrix,float *inversematrix,int n)
{
    int det=determinant(matrix,n);
    if(det==0)
    {
        return 0;
    }
    int *cofactormatrix=(int *)malloc(n*n*sizeof(int));
    int *adjointmatrix=(int *)malloc(n*n*sizeof(int));
    
    findcofactor(matrix,cofactormatrix,n);
    findadjoint(cofactormatrix,adjointmatrix,n);
    
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            inversematrix[i*n+j]=(float)adjointmatrix[i*n+j]/det;
        }
    }
    
    free(cofactormatrix);
    free(adjointmatrix);
    
    return 1;
}

int main()
{
    int m,n,p,q,k;
    char str[10];

    // file for printing result
    FILE *file3 =fopen("result.txt","w");
    if(file3==NULL) 
    {
        printf("Error in opening file.txt\n");
        return 1;
    }

    printf("Enter the type of operation: ");
    scanf("%s",str);  

    if(strcmp(str,"one")==0)  // string comparison
    {
        printf("Enter dimensions of matrix(rows cols): ");
        scanf("%d%d",&m,&n);

        // for single operation allocation of array size
        int *A=(int *)malloc(MAX*MAX*sizeof(int));

        FILE *file1=fopen("matrix1.txt","w");
        if(file1==NULL) 
        {
            printf("Error opening the file.txt\n");
            free(A);
            return 1;
        }

        printf("Enter values for matrix:\n");
        for(int i=0; i< m;i++) 
        {
            for(int j=0;j<n;j++) 
            {
                scanf("%d",&A[i* n+j]); // taking the input
                fprintf(file1,"%d ",A[i*n+j]); // printing the elements
            }
            fprintf(file1,"\n");
        }
        fclose(file1);

        getchar();  //newline
        char ch1;
        printf("Enter the operation (D,T,t,C,A): ");
        scanf("%c",&ch1); 

        switch(ch1) 
        {
        case 'D':
            if (m!=n) 
            {
                printf("The number of rows and columns should be equal\n");
            } 
            else 
            {
                int det=determinant(A,n);
                fprintf(file3,"Determinant:%d\n",det);  // printing the value
            }
            break;
        case 'T':
        {
            int transpose1[m][n]; // for transpose
            transposingthematrix((int *)A,(int *)transpose1,m,n);
            for(int i=0;i<n;i++) 
            {
                for(int j=0;j<m;j++) 
                {
                    fprintf(file3,"%d ",transpose1[i][j]);  // printing the matrix
                }
                fprintf(file3,"\n");
            }
        }
        break;
        case 't':
        {
            traceofmatrix(A,n,file3);
        }
            break;
        case 'C':
        {
            int *cofactormatrix=(int *)malloc(n*n*sizeof(int));
            findcofactor(A,cofactormatrix,n);
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<n;j++)
                {
                    fprintf(file3,"%d",cofactormatrix[i*n+j]);
                }
                fprintf(file3,"\n");
            }
            free(cofactormatrix);
        }
            break;
        case 'A':
        {
            int *cofactormatrix=(int *)malloc(n*n*sizeof(int));
            findcofactor(A,cofactormatrix,n);
            
            int *adjointmatrix=(int *)malloc(n*n*sizeof(int));
            findadjoint(cofactormatrix,adjointmatrix,n);
            
            fprintf(file3,"Adjoint matrix:\n");
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<n;j++)
                {
                fprintf(file3,"%d",adjointmatrix[i*n+j]);
                }
                fprintf(file3,"\n");
            }
        }
            break;
        case 'I':
        {
            float *inversematrix=(float *)malloc(m*n*sizeof(float));
            int result=findinverse(A,inversematrix,n);
            if(result)
            {
                fprintf(file3,"Inversematrix:\n");
                for(int i=0;i<n;i++)
                {
                    for(int j=0;j<n;j++)
                    {
                        fprintf(file3,"%.2f ",inversematrix[i*n+j]);
                    }
                }
                fprintf(file3,"\n");
            }
            else
            {
                fprintf(file3,"inverse is invalid");
            }
        }
        break;
        default:
        {
            printf("Invalid operation.\n");
        }
    }
    }
    else if(strcmp(str,"two")==0)  // string comparison
    {
        printf("Enter dimensions of matrix 1 (rows cols): ");
        scanf("%d%d",&m,&n);
        printf("Enter dimensions of matrix 2 (rows cols): ");
        scanf("%d%d",&p,&q);

        int *A = (int *)malloc(m*n*sizeof(int));
        int *B = (int *)malloc(p*q*sizeof(int));

        int *result1 = (int *)malloc(m*n*sizeof(int)); // matrix size allocation for '+' and '-'
        int *result2 = (int *)malloc(m*q*sizeof(int)); // matrix size allocation for '*'

        if (A==NULL||B==NULL||result1==NULL||result2==NULL) 
        {
            printf("Memory allocation failed.\n");
            return 1;
        }

        FILE *file1 =fopen("matrix1.txt","w");
        if (file1==NULL) 
        {
            printf("Error opening matrix1.txt\n");
            free(A);free(B);free(result1);free(result2);
            return 1;
        }

        FILE *file2 = fopen("matrix2.txt","w");
        if (file2==NULL) 
        {
            printf("Error opening matrix2.txt\n");
            free(A);free(B);free(result1);free(result2);
            fclose(file1);
            return 1;
        }

        printf("Enter values for matrix 1:\n");
        for(int i=0;i<m;i++) 
        {
            for (int j=0;j<n;j++) 
            {
                scanf("%d",&A[i*n+j]); // taking the input
                fprintf(file1,"%d ",A[i*n+j]); // printing the matrix1
            }
            fprintf(file1,"\n");
        }
        fclose(file1);

        printf("Enter values for matrix 2:\n");
        for (int i=0;i<p;i++) 
        {
            for (int j=0;j<q;j++) 
            {
                scanf("%d",&B[i * q + j]); // taking the input
                fprintf(file2,"%d ",B[i * q + j]); // printing the matrix2
            }
            fprintf(file2,"\n");
        }
        fclose(file2);

        getchar();  // To clear the buffer
        char ch2;
        printf("Enter the operation(+, -, *): ");
        scanf("%c",&ch2);  // Fix for getchar issue

        switch(ch2) 
        {
            case '+':
                addition(A,B,result1,m,n,p,q,file3);
                break;

            case '-':
                subtraction(A,B,result1,m,n,p,q,file3);
                break;

            case '*':
                multiplication(A,B,result2,m,n,p,q,file3);
                break;
            
            default:
                printf("Invalid operation.\n");
        }

        free(A);
        free(B);
        free(result1);
        free(result2);
    }
    else
    {
        printf("Invalid operation type.\n");
    }

    fclose(file3);
    return 0;
}
