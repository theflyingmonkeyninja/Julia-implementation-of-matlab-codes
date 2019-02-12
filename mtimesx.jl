function mtimesx(A,B,moperandi)

    ## function to compute the matrix multiplication for multidimensional matrices
    ## A and B are two matrices have same dimesnion
    ## moperandi is the variable that defines the how the operation will be performed
    ## moperandi is a string and could be either ""vectorised" indicating use of simd
    ## or "threaded" indicating use of multithread approach the default way the function
    ## runs as single threaded operation


    fmatshap = collect(size(A));
    fmatshap = fmatshap[3:end];
    fmatshap = reshape(fmatshap,1,length(fmatshap));



    if ndims(A)==ndims(B)
        k=1;
        for i=3:ndims(A)
            if size(A,i) != size(B,i)
                error("Matrices/vectors are incorrectly defined")
            end

        end
    else
        error("Matrices/vectors are incorrectly defined")
    end
        ## getting the size of the matrices
    sa1 = size(A,1);
    sa2 = size(A,2);
    sb1 = size(B,1);
    sb2 = size(B,2);

    ## converting to 3d matrices
    A = reshape(A,sa1,sa2,prod(fmatshap));
    B = reshape(B,sb1,sb2,prod(fmatshap));

    ## getting the size of the matrices
    (sa1,sa2,sa3) = size(A);
    (sb1,sb2,sb3) = size(B);

    ## computing of the matrix
    if sa3 == sb3 && sb1 == sa2

        if moperandi == "vectorised" ## using simd
             C = zeros(sa1,sb2,sa3);

            @inbounds @simd for i = 1:sb3
                C[:,:,i] = A[:,:,i]*B[:,:,i];
            end

        elseif moperandi == "threaded" ## using multithreading
            C = zeros(sa1,sb2,sa3);

            Threads.@threads for i = 1:sb3
                C[:,:,i] = A[:,:,i]*B[:,:,i];
            end

        else
            C = zeros(sa1,sb2,sa3);

            for i = 1:sb3 ## single threaded
                C[:,:,i] = A[:,:,i]*B[:,:,i];
            end
        end
    else

        error("Matrices/vectors are incorrectly defined")

    end

    ## reconverting the matrix to the correct shape
    spvect = zeros(Int64,1,2+length(fmatshap));
    spvect[1] = sa1; spvect[2] = sb2;

    for i = 1:length(fmatshap)
        spvect[2+i] = fmatshap[i];
    end


    C = reshape(C,Tuple(spvect));# reshape takes tuples as second argument
    return C;


end
