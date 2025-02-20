module intro 
        implicit none

        !for-loops variables.
        integer*4 ii, jj, kk

        !dimensions of the matrices. 
        integer*4 n_rows_AA, n_rows_BB, n_columns_AA, n_columns_BB !numbers of rows and columns of the 2 matrices
        
        !cpu_time variables.
        real*4 start, finish
        real*4 max_error_allowed

        !matrices
        real*4, allocatable, dimension(:,:) :: AA, BB, CC_1, CC_2, CC_matmul

    end module intro

    !matrix multiplication

    program Ex3 
    
        use intro
        implicit none
        
    !enter the matrices dimension (by keyboard).
        READ(*,*) n_rows_AA, n_columns_AA
        READ(*,*) n_rows_BB, n_columns_BB
        ALLOCATE(AA(n_rows_AA, n_columns_AA), BB(n_rows_BB, n_columns_BB))

    !check if the matrices are compatibles for the product.
        IF (n_columns_AA /= n_rows_BB) THEN
            PRINT*, 'In order to do the matrix product the #columns of A has to be equal to #rows of B'
            return
        END IF
       
        
    !CHOICE OF FILLING THE MATRIX ENTRIES:

        
    ! Random filling:

        CALL random_number(AA) !type output should depend on the variable's type's input
        CALL random_number(BB)
        

    !MATRIX MULTIPLICATION IN THREE WAYS (with cpu_time measurement):
    
    !1) 3 for-loops, Usual Matrix Product: 
        
        ALLOCATE(CC_1(n_rows_AA, n_columns_BB))
        CC_1 = 0.0

            
        CALL cpu_time(start)

            DO ii = 1, n_rows_AA
                DO jj = 1, n_columns_BB
                    DO kk = 1, n_columns_AA 
                        CC_1(ii,jj) = CC_1(ii,jj) + AA(ii, kk) * BB(kk, jj)
                    ENDDO
                ENDDO
            ENDDO

        CALL cpu_time(finish)

        PRINT*, finish - start

            
    !2) 3 for-loops with inverted indices:
        
        ALLOCATE(CC_2(n_rows_AA, n_columns_BB))
        CC_2 = 0.0
            

            
        CALL cpu_time(start)

        !3 for-loops algorithm with inverted indeces
            DO kk = 1, n_columns_AA
                DO ii = 1, n_rows_AA
                    DO jj = 1, n_columns_BB
                        CC_2(ii,jj) = CC_2(ii,jj) + AA(ii, kk) * BB(kk, jj)
                    ENDDO
                ENDDO
            ENDDO

        CALL cpu_time(finish)

        PRINT*, finish - start

    !3) MATMUL product:

        ALLOCATE(CC_matmul(n_rows_AA, n_columns_BB))
        CC_matmul = 0.0


        CALL cpu_time(start)
        CC_matmul = matmul(AA,BB)
        CALL cpu_time(finish)

        PRINT*, finish - start


       
    !Deallocate matrices
        DEALLOCATE(CC_1, CC_2)
        DEALLOCATE(AA, BB, CC_matmul)

    end program Ex3

    subroutine print_matrix(AA,nn,mm)
        integer*4 nn, mm, kk
        real*4 AA(nn,mm) !nn = #rows, mm = #columns
        DO kk=1,nn
        PRINT '(20f6.2)', AA(kk,1:mm)
        ENDDO
    endsubroutine   

   subroutine comparing_matrix(AA, BB, nn, mm, eps)
        integer*4 nn, mm, ii, jj, err_counts
        real*4 AA(nn, mm), BB(nn, mm), eps
        err_counts = 0

        DO ii = 1, nn
            DO jj = 1, mm
                IF(ABS(AA(ii,jj) - BB(ii,jj)) > eps) THEN
                    PRINT*, 'Entry (',ii, ',', jj,')', 'of the matrices are different'
                    err_counts = err_counts + 1
                ENDIF
            ENDDO
        ENDDO

        IF(err_counts /= 0) THEN
            PRINT*, new_line('a'), 'There are', err_counts, 'entries in the two matrices that differ more than a value of', eps
        ELSE
            PRINT*, 'Matrices result the same up to a value', eps
        ENDIF

    endsubroutine 