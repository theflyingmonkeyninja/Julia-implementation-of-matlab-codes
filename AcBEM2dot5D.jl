## Inputs
# Ub displacements matrix of size (NW,Nkx,Nnodes);
# kx wavenumbers matrix of (nw,nkx);
# wm moving frequency
# ws static frequency matrix of (nw,nkx)
# route needs to be string specifying calculation files

## Outputs
# Pf is pressure on field points matrix of dimension (Nevals,Nkx,Nw)
# Pb is pressure on boundary points matrix of dimension (Nnodes,Nkx,Nw)

## OPts
# BEdomain type of domain needs to be 'interior' for interior and [] for exterior
# sampling needs to be moving for moving sampling and [] for static
# c is wavespeed c == [] implies c = 340 m/s without damping


## The solver function

function AcBEMsolver(route,c,n,w,kx,Ub,FPnt,boundarypressure,BEdomain,sampling,Opts...)
    ## The boundary Element solver


    if isempty(route)==1 # testting purpose

        ConM = Opts[1];  Nodes = Opts[2];

    else

        # reading file
        (~,conect_mat,~,~,~,~,mat_IntBEM,BGNodes)=readcalcfile(route);
        # making connectivity matrix
        (IntGConM) = Connectivity(mat_IntBEM,BGNodes,conect_mat);
        # making appropriate connectivity
        (ConM,~,Nodes) = dictionary(IntGConM,BGNodes);

    end
    ## Initilising
    (xi,wf,phi,Y,Z,Ny,Nz,jacobi,nodes,Vs,Vs0,LPnt) = AcBEMInit(c,n,ConM,Nodes,BEdomain);

    if size(FPnt,2)>2
        FPnt = FPnt[:,2:3];
    end


    Nw = size(w,1); Nkx = size(kx,2);
    Neval = size(FPnt,1); Nnodes = size(Nodes,1);

    Pf = zeros(ComplexF64,Neval,Nkx,Nw);

    if boundarypressure == "yes"

        Pb = zeros(ComplexF64,Nnodes,Nkx,Nw);

    else

        Pb = [];

    end

  for nw = 1:Nw


        if sampling=="moving"
            # computing boundary matrices for moving source
            (matHb,matGb,matHf,matGf) = AcBEBndMatrices(c,n,w[nw,:],kx[nw,:],xi,wf,phi,Y,Z,Ny,Nz,Vs0,jacobi,nodes,LPnt,FPnt,ConM,Nodes,Vs);

            Vb = permutedims(repeat(permutedims(w[nw,:].^2,3,2,1),inner=(size(Ub,1)),1,1).*Ub(nw,:,:),3,2,1);

        else
            # computing boundary matrices for non-moving source
            (matHb,matGb,matHf,matGf) = AcBEBndMatrices(c,n,w[nw],kx[nw,:],xi,wf,phi,Y,Z,Ny,Nz,Vs0,jacobi,nodes,LPnt,FPnt,ConM,Nodes,Vs);

            Vb = w[nw].*Ub;
        end

        rho = 1.24;
        Vb = (1im*1im)*rho*Vb;

        if boundarypressure=="yes"
            bpressure=1;
        else
            bpressure=0;
        end

        # computing solution

        for nkx = 1:Nkx

            if bpressure==1

                (pf,pb) = AcBEMsolution(matHf,matGf,matGb,matHb,nkx,Vb);

                Pf[:,nkx,nw] = pf;
                Pb[:,nkx,nw] = pb;

            else
                (pf,~) = AcBEMsolution(matHf,matGf,matGb,matHb,nkx,Vb);

                Pf[:,nkx,nw] = pf;
            end

        end
    end

    return Pf,Pb;
end

## Initilising BEM
function  AcBEMInit(c,n,ConM,Nodes,BEdomain)
    ## function to inilise BEM

    # Getting geometry variables as a function of xi
    (xi,wf,phi,Y,Z,Ny,Nz,jacobi,LPnt,nodes) = geometryvariables(n,ConM,Nodes);

    # Locating singularity
    (~,ElemId,~) = singLocate(LPnt,ConM,Nodes,nodes);

    # Daigonal Teradm for acoustics
    (Hs) = AcDiagTeradm(c,n,wf,Y,Z,Ny,Nz,jacobi,LPnt,ConM,ElemId,BEdomain);

    # Relative distances for static tractions
    (~,ys,zs,~,nys,nzs) = singDist([],Nodes,ConM,ElemId,LPnt,xi);

    # Static Tractions
    (~,Hs0) = AcFullSpace2dot5D(c,0,ys,zs,0,nys,nzs);

    return xi,wf,phi,Y,Z,Ny,Nz,jacobi,nodes,Hs,Hs0,LPnt;
end

function AcDiagTeradm(c,n,wf,Y,Z,Ny,Nz,jacobi,LPnt,LConM,ElemId,BEdomain)
    ## function for computing daigonal teradm for acoustics

    wf = repeat(wf',inner=(1,1,size(LConM,1)));
    He = zeros(size(LPnt,1),1,size(LConM,1));


    (yrel,zrel,~,rrel,ny,nz) = relativedistance([],LPnt,Y,Z,Ny,Nz,jacobi);
    (~,H) = AcFullSpace2dot5D(c,0,yrel,zrel,0,ny,nz); #computing greens function for all relative distance

    for i = 1:size(LPnt,1)

        indx=collect(size(LConM,1)*n*(i-1).+(collect(1:size(LConM,1)*n)));

        radm = rrel[indx];

        (~,ibs) = ismembertol(radm,rrel,1e-6);

        Hs4D = reshape(H,1,1,1,length(rrel));
        H0 = reshape(Hs4D[:,:,:,ibs],1,1,size(LConM,1)*n);
        H0 = real(H0);   # real part of static function

        H0r = reshape(H0[1,1,:],n,1,size(LConM,1)) ;
        Hre = zeros(1,1,size(LConM,1)) ;

        Hre[:,:,:] = mtimesx(permutedims(wf,[2,1,3]),H0r[:,:,:].*jacobi[:,:,:],"") ;

        for j=1:size(ElemId,1)

            Hre[1,1,ElemId[j,i]] = 0;

        end

        He[i,:,:] = Hre;

    end

    Hs = reshape(sum(He,dims=3),size(LConM,1),1);

    if BEdomain !="interior"
        Idn = ones(size(LConM,1),1);
        Hs = Idn-Hs;
    else
        Hs = -Hs;
    end

    return Hs;
end

## BEM Solution

function AcBEMsolution(matHf,matGf,matGb,matHb,Nkx,Vb)
    ## function to compute pressure at fieldpoint

    MatHb = matHb[:,:,Nkx]; MatGb = matGb[:,:,Nkx];
    MatHf = matHf[:,:,Nkx]; MatGf = matGf[:,:,Nkx];
    Vb = Vb[:,Nkx];

    # Getting pressure on the the boundary
    Pb = MatHb\MatGb*Vb;

    # Getting the displacement at the field Points
    Pf = -(MatHf*Pb-MatGf*Vb);

    return  Pf,Pb;

end

function AcIntRepEqn(Pn,Vn,LConM,LNodes,FPnt,nkx)
    ## function to assemble the Green's Traction and Displacement matrices for fieldpoint

    MatG = zeros(ComplexF64,size(FPnt,1),size(LNodes,1),1,nkx);
    MatH = zeros(ComplexF64,size(FPnt,1),size(LNodes,1),1,nkx);

    Pn = reshape(Pn,size(FPnt,1),2,size(LConM,1),nkx);
    Vn = reshape(Vn,size(FPnt,1),2,size(LConM,1),nkx);


    # Assembling the Green's Traction and Displacement Matrices
    for i = 1:size(LConM,1)

        MatG[:,[LConM[i,2], LConM[i,3]],:,:] = MatG[:,[LConM[i,2], LConM[i,3]],:,:]+Pn[:,:,i,:];
        MatH[:,[LConM[i,2], LConM[i,3]],:,:] = MatH[:,[LConM[i,2], LConM[i,3]],:,:]+Vn[:,:,i,:];

    end

    MatG=MatG[:,:,:];
    MatH=MatH[:,:,:];

    return MatG,MatH;
end

function AcIntRepEqnMat(c,n,w,kx,wf,phi,Y,Z,Ny,Nz,jacobi,FPnt,LConM)

    ## function to make Green's Traction and Dispalacement  matrices for fieldpoint

    nkx = length(kx); kx= kx'; w=w'; # f and kx must be a column vector
    Jacobi = jacobi; wf = repeat(wf,inner=(1,1,size(LConM,1),nkx)); phi =phi;

    Vn = zeros(ComplexF64,size(FPnt,1),2,size(LConM,1)*nkx);
    Pn = zeros(ComplexF64,size(FPnt,1),2,size(LConM,1)*nkx);

    (yrel,zrel,~,rrel,ny,nz) = relativedistance([],FPnt,Y,Z,Ny,Nz,jacobi);
    (Ps,Us,rs) = AcFullSpace2dot5DBEM(c,kx,yrel,zrel,w);

     for i = 1:size(FPnt,1)

        indx=collect(size(LConM,1)*n*(i-1).+(collect(1:size(LConM,1)*n)));

        radm = rrel[indx];

        (~,ibs) = ismembertol(radm,rs,1e-6);

        Ps4D = reshape(Ps,1,1,nkx,length(rs));
        Vs4D = reshape(Us,1,1,nkx,length(rs));

        P = reshape(Ps4D[:,:,:,ibs],1,1,nkx*size(LConM,1)*n);
        U = AcFullSpace2dot5DBEMvelocity(Vs4D[:,:,:,ibs],yrel,zrel,radm,ny,nz);


        Pr = reshape(P[1,1,:],nkx,n,1,size(LConM,1));
        Vr = reshape(U[1,1,:],nkx,n,1,size(LConM,1));

        Pr = permutedims(Pr,[2 3 4 1]);
        Vr = permutedims(Vr,[2 3 4 1]);

        # vectors for greens displacement on elements
        Pr1 = zeros(ComplexF64,1,1,size(LConM,1),nkx);
        Pr2 = zeros(ComplexF64,1,1,size(LConM,1),nkx);

        # vectors for greens tractions on nodes of element
        Vr1 = zeros(ComplexF64,1,1,size(LConM,1),nkx);
        Vr2 = zeros(ComplexF64,1,1,size(LConM,1),nkx);

        # Greens displacemnet matrix element computed at first nodal point for all elements
        Pr1[:,:,:,:] = mtimesx((wf),(phi[:,1].*Pr.*Jacobi),"");
        # Greens displacemnet matrix element computed at second nodal point for all elements
        Pr2[:,:,:,:] = mtimesx((wf),(phi[:,2].*Pr.*Jacobi),"");

        # Greens traction matrix element computed at first nodal point for all elements
        Vr1[:,:,:,:] = mtimesx((wf),(phi[:,1].*Vr.*Jacobi),"");
        # Greens traction matrix element computed at second nodal point for all elements
        Vr2[:,:,:,:] = mtimesx((wf),(phi[:,2].*Vr.*Jacobi),"");

        # converting displacement in 3d matrix for speed
        Pr1=Pr1[:,:,:];
        Pr2=Pr2[:,:,:];

        # converting traction in 3d matrix for speed
        Vr1=Vr1[:,:,:];
        Vr2=Vr2[:,:,:];

        # the final matrices of pressure and displacement
        Pn[i,:,:] = [Pr1  Pr2];

        Vn[i,:,:] = [Vr1  Vr2];
    end

    return  Pn, Vn;
end

## Boundary Solution Matrices

function AcBEBndMatrices(c,n,w,kx,xi,wf,phi,Y,Z,Ny,Nz,Vs0,jacobi,nodes,LPnt,FPnt,LConM,LNodes,Vs)

    ## function to make BE matrices

    # Elements for BE matrices
    (Vn,Pn,CPntId,Lpnt) = AcBndIntEqnMat(c,n,w,kx,xi,wf,phi,Y,Z,Ny,Nz,Vs0,jacobi,nodes,LPnt,LConM,LNodes);

    #Size of kx vector
    nkx = length(kx);

    # Boundary Matrices for Greens traction and displacement on Boundary
    (MatGb,MatHb) = AcBndIntEqnInt(Vs,Pn,Vn,Lpnt,CPntId,LNodes,LConM,nkx);

    if !isempty(FPnt)
        # Getting Greens traction and displacement for fieldpoint
        (Pn,Vn)=AcIntRepEqnMat(c,n,w,kx,wf,phi,Y,Z,Ny,Nz,jacobi,FPnt,LConM);

        # Building the traction and displacement matrices for fieldpoint
        (MatGf,MatHf)=AcIntRepEqn(Pn,Vn,LConM,LNodes,FPnt,nkx);

    else

        MatHf=[]; MatGf=[];

    end

    return MatHb,MatGb,MatHf,MatGf;
end

function AcBndIntEqnInt(Vs,Pn,Vn,Lchf,CPntId,LNodes,LConM,nkx)
    ## function to assemble the Greens Traction and Displacement matrices on boundray

    MatG = zeros(ComplexF64,Lchf,size(LNodes,1),1,nkx);
    MatH = zeros(ComplexF64,Lchf,size(LNodes,1),1,nkx);


    # Assembling the Green's Traction and Displacement Matrices
    MatI = zeros(Lchf,size(LNodes,1));

    for i = 1:size(LConM,1)

        MatG[:,[LConM[i,2], LConM[i,3]],:,:] = MatG[:,[LConM[i,2], LConM[i,3]],:,:]+Pn[:,:,i,:];
        MatH[:,[LConM[i,2], LConM[i,3]],:,:] = MatH[:,[LConM[i,2], LConM[i,3]],:,:]+Vn[:,:,i,:];

    end

    for i = 1:size(LNodes,1)

        MatI[i,CPntId[i]] = MatI[i,CPntId[i]]+Vs[i];

    end

    MatH=MatH[:,:,:];
    MatG=MatG[:,:,:];
    MatH = MatH+repeat(MatI,inner=(1,1,nkx));

    return MatG,MatH;
end


function AcBndIntEqnMat(c,n,w,kx,xi,wf,phi,Y,Z,Ny,Nz,Us0,jacobi,nodes,LPnt,LConM,LNodes)
    ## function to make the Greens Traction and Displacement matrices on boundary

    nkx = length(kx); kx = kx'; w=w'; # f and kx needs to be a column vector
    Jacobi=jacobi;  wf = repeat(wf,inner=(1,1,size(LConM,1),nkx)); phi = phi;

    (nodeID,ElemId,CPntId) = singLocate(LPnt,LConM,LNodes,nodes);
    l = size(LPnt,1);

    Lpnt = size(LPnt,1);

    Pn = zeros(ComplexF64,Lpnt,2,size(LConM,1)*nkx);
    Vn = zeros(ComplexF64,Lpnt,2,size(LConM,1)*nkx);

    (yrel,zrel,~,rrel,ny,nz) = relativedistance([],LPnt,Y,Z,Ny,Nz,jacobi);
    (G,H,rs) = AcFullSpace2dot5DBEM(c,kx,yrel,zrel,w);

    for i = 1:size(LPnt,1)

        indx=collect(size(LConM,1)*n*(i-1).+(collect(1:size(LConM,1)*n)));

        radm = rrel[indx];

        yrels = yrel[indx];
        zrels = zrel[indx];

        nys = ny[indx];
        nzs = nz[indx];

        (~,ibs) = ismembertol(radm,rs,1e-6);

        Ps4D = reshape(G,1, 1, nkx, length(rs));
        Vs4D = reshape(H,1, 1, nkx, length(rs));

        P = reshape(Ps4D[:,:,:,ibs],1, 1, nkx*size(LConM,1)*n)
        V = AcFullSpace2dot5DBEMvelocity(Vs4D[:,:,:,ibs],yrels,zrels,radm,nys,nzs);

        # reshaping greens displacement vector
        Pr = reshape(P[1,1,:],nkx,n,1,size(LConM,1)) ;
        Pr = permutedims(Pr,[2 3 4 1]);

        # reshaping greens traction vector
        Vr = reshape(V[1,1,:],nkx,n,1,size(LConM,1)) ;
        Vr = permutedims(Vr,[2 3 4 1]);

        # vectors for greens displacements on nodes of element
        Pr1 = zeros(ComplexF64,1,1,size(LConM,1),nkx) ;
        Pr2 = zeros(ComplexF64,1,1,size(LConM,1),nkx) ;

        # vectors for greens tractions on nodes of element
        Vr1 = zeros(ComplexF64,1,1,size(LConM,1),nkx) ;
        Vr2 = zeros(ComplexF64,1,1,size(LConM,1),nkx) ;

        ## Regular Integration

        # Greens displacemnet matrix element computed at first nodal point for all elements
        Pr1[:,:,:,:] = mtimesx((wf),(phi[:,1].*Pr.*Jacobi),"");
        # Greens displacemnet matrix element computed at second nodal point for all elements
        Pr2[:,:,:,:] = mtimesx((wf),(phi[:,2].*Pr.*Jacobi),"");

        # Greens traction matrix element computed at first nodal point for all elements
        Vr1[:,:,:,:] = mtimesx((wf),(phi[:,1].*Vr.*Jacobi),"");
        # Greens traction matrix element computed at second nodal point for all elements
        Vr2[:,:,:,:] = mtimesx((wf),(phi[:,2].*Vr.*Jacobi),"");

        ## singular Integration
        if i<l+1

            for j=1:size(ElemId,1)

                (mphi) = modifiedshapefunction(xi);
                mphi = repeat(mphi,1,1,1,nkx);

                k=2*(i-1)+j;

                V0ss = Us0[:,:,collect((k-1)*n+1:k*n)];
                Vss = real(V0ss); # real part for static function
                Vms = real(V0ss); # real part for modified static function

                Vssr = reshape(Vss[1,1,:],n,1,1) ;
                Vssr = repeat(Vssr,1,1,1,nkx);

                Vmsr = reshape(Vms[1,1,:],n,1,1) ;
                Vmsr = repeat(Vmsr,1,1,1,nkx);

                Vdr1 = mtimesx((wf[:,:,ElemId[j,i],:]),(phi[:,1].*(Vr[:,:,ElemId[j,i],:]-Vssr).*Jacobi[:,:,ElemId[j,i],:]),"") ;
                Vdr2 = mtimesx((wf[:,:,ElemId[j,i],:]),(phi[:,2].*(Vr[:,:,ElemId[j,i],:]-Vssr).*Jacobi[:,:,ElemId[j,i],:]),"") ;

                V0sr1 = mtimesx((wf),(mphi[:,1].*(Vmsr).*Jacobi),"") ;
                V0sr2 = mtimesx((wf),(mphi[:,2].*(Vmsr).*Jacobi),"") ;

                # singular Greens tractions
                #Vr1[1,1,ElemId[j,i],:] = Vdr1+V0sr1 ;
                #Vr2[1,1,ElemId[j,i],:] = Vdr2+V0sr2 ;
            end
        else
            break;
        end
        ##  Assembling the elements

        Pr1=Pr1[:,:,:];
        Pr2=Pr2[:,:,:];

        Vr1=Vr1[:,:,:];
        Vr2=Vr2[:,:,:];

        Pn[i,:,:] = [Pr1  Pr2];

        Vn[i,:,:] = [Vr1  Vr2 ];

    end

    return Vn,Pn,CPntId,Lpnt;
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------
## BEM functions
#-------------------------------------------------------------------------------------------------------------------------------------------------------
function singLocate(LPnt,LConM,LNodes,nodes)
    ## Function that loacates the singularity for both fullspace and halfspace

    CPntId = zeros(Int8,1,size(LNodes,1));
    NdeId = zeros(Int8,2,size(LNodes,1)); # node id of singularity
    ElemId = zeros(Int8,2,size(LConM,1)); # element id of singularity

    for i = 1: size(LPnt,1)

        # finding singular element and nodal location of singularity on that element
        Pxy = ones(2,1)*transpose(LPnt[i,:]);
        pntmat = repeat(Pxy,inner=[1,1,size(LConM,1)]);
        testpnt = sum(abs.(nodes[:,:,:]-pntmat),dims=2);
        (nodeid,elmtid) = findvals(testpnt,0);

        NdeId[:,i] = nodeid;  # node id of load points
        ElemId[:,i] = elmtid; # element id of load point

        CPntId[i] = LConM[elmtid[1],nodeid[1]+1]; # location of singularity
    end
    return NdeId,ElemId,CPntId;
end

function singDist(Method,LNodes,LConM,ElemId,LPnt,xi)
    ## Function to get the relative distances on singular element for fullspace solution

    zf = zeros(size(LPnt,1),1); # source point cordinate
    ys = zeros(2*length(xi),size(LPnt,1)); # singular distances in y
    zs = zeros(2*length(xi),size(LPnt,1)); # singular distances in z

    nys = zeros(2*length(xi),size(LPnt,1)); # noradmals on singular element
    nzs = zeros(2*length(xi),size(LPnt,1)); # noradmals on singular element

    for i= 1:size(LPnt,1)

        y_sing = zeros(2,length(xi));
        z_sing = zeros(2,length(xi));

        ny_sing = zeros(2,length(xi));
        nz_sing = zeros(2,length(xi));

        for j=1:size(ElemId,1)

            NCordinates = LNodes[LConM[ElemId[j,i],2:3],2:end];
            (y,z,ny,nz) = isoparametric(NCordinates,xi);

            y = y.-LPnt[i,1];

            if Method!="hlfspc"
                z = z.-LPnt[i,2];
            end

            y_sing[j,:] = y;
            z_sing[j,:] = z;

            ny_sing[j,:] = ny;
            nz_sing[j,:] = nz;

        end

        y_sing = y_sing';
        z_sing = z_sing';

        ny_sing = ny_sing';
        nz_sing = nz_sing';

        ys[:,i] = y_sing[:];
        zs[:,i] = z_sing[:];

        nys[:,i] = ny_sing[:];
        nzs[:,i] = nz_sing[:];
    end

    ys = ys[:];
    zs = zs[:];

    nys = nys[:];
    nzs = nzs[:];

    rs = sqrt.(ys.^2+zs.^2);

    if Method=="hlfpsc'"  # for halfspace cordinate of the force points
        zf = LPnt[:,2];
    end

    return zf,ys,zs,rs,nys,nzs;
end

function relativedistance(Method,LPnt,Y,Z,Ny,Nz,jacobi)
    ## Function that computes the relative distance required for BEM matrices

    y = zeros(size(Y[:],1),size(LPnt,1));
    z = zeros(size(Z[:],1),size(LPnt,1));


    for i=1:size(LPnt,1)

        y[:,i] = Y[:].-LPnt[i,1];

        if Method =="hlfspc"

            z[:,i] = Z[:];

        else

            z[:,i] = Z[:].-LPnt[i,2];
        end
    end

    ny = Ny./jacobi; ny = ny[:];
    nz = Nz./jacobi; nz = nz[:];

    ny = repeat(ny,inner=(size(LPnt,1),1));
    nz = repeat(nz,inner=(size(LPnt,1),1));

    mat =[y[:] z[:]];

    YZrel = round.(mat,digits=5);

    #YZrel = unique(YZrel,'rows');

    Yrel = YZrel[:,1]; # distance in y
    Zrel = YZrel[:,2]; # distance in z

    rrel= sqrt.(Yrel.^2+Zrel.^2);

    if Method=="hlfspc"
        Zf = LPnt[:,2];
    else
        Zf = [];
    end
    return Yrel,Zrel,Zf,rrel,ny,nz;
end

#-------------------------------------------------------------------------------------------------------------------------------------------------------
## Geometry Functions BEM
#-------------------------------------------------------------------------------------------------------------------------------------------------------

function geometryvariables(n,LConM,LNodes)
    ## function to get geometry variables in xi cordinate system

    # integration point and weights
    (xi,wf) = gaussrule(n);

    # shapefunction for isoparametric discreetisation
    (phi,~) = shapefunction(xi);

    # geometry variables as function of xi
    (X,Z,Nx,Nz,jacobi,LPnt,nodes) = geometrytransforadm(LConM,LNodes,xi);

    return xi,wf,phi,X,Z,Nx,Nz,jacobi,LPnt,nodes;
end

function geometrytransforadm(LConM,LNodes,xi)
    ## Returns Y,Z;Ny,Nz; jacobi as a function of xi for all elements
    # Y is a 3d matrix of [N,1,size(elements,1)] containg y(xi) for all elements
    # Z is a 3d matrix of [N,1,size(elements,1)] containg z(xi) for all elements
    # NY is a 3d matrix of [N,1,size(elements,1)] containg ny(xi) for all elements
    # Nz is a 3d matrix of [N,1,size(elements,1)] containg nz(xi) for all elements
    # Jacobi is a 3d matrix of [N,1,size(elements,1)] containg jacobi(xi) for all elements

    # initialising matrix
    NCords = zeros(size(LConM,2)-1,size(LNodes,2)-1,size(LConM,1));

    Y = zeros(length(xi),1,size(LConM,1));
    Z = zeros(length(xi),1,size(LConM,1));

    Ny = zeros(length(xi),1,size(LConM,1));
    Nz = zeros(length(xi),1,size(LConM,1));

    # getting the cordinates associated to nodesID of an element and storing them in 3d matrix
    # node_cordinate: 1st row has the cordinates of 1st nodeID of element 1  while second row has cordinates of 2nd nodeID of element1
    # node_cordinate: 3d dimension does the same for all elements


    for j = 1:size(LConM,1)
        NCords[:,:,j] = LNodes[LConM[j,2:3],2:end];
    end

    LPnt = LNodes[:,2:end];
    Nodes = NCords;

    for j = 1:size(NCords,3)

        (y,z,ny,nz) = isoparametric(NCords[:,:,j],xi);

        Y[:,:,j] = y;
        Z[:,:,j] = z;

        Ny[:,:,j] = ny;
        Nz[:,:,j] = nz;
    end

    jacobi = sqrt.(Ny.^2+Nz.^2);

    return Y,Z,Ny,Nz,jacobi,LPnt,Nodes;
end

function shapefunction(xi)

    ## returns the shapefunction and its derivate for linear element
    # xi is a vector of size [n,1]
    # phi,dphi are matrix of size [n,2]

    phi = [-0.5*(xi.-1) 0.5*(xi.+1)];
    dphi = [-0.5.*ones(size(xi,1)) 0.5.*ones(size(xi,1))];
    return phi,dphi;
end

function isoparametric(nodes_cordinates,xi)

    ## Returns cordinates y,z; shape function phi; noradmals nx,ny; and jacobi as a function of xi
    # nodes_cordinates is matrix of size [2,2] the columns of which represent cordinates x,y
    # xi is vector of size[n,1]
    # y,z,ny,nz,jacobi are vectors of size[n,1]
    # phi is matrix of size[n,2]

    # Linear Shape functions and their derivatives
    (phi,dphi) = shapefunction(xi);

    # Local cordinate of x and y
    y = phi*nodes_cordinates[:,1];
    z = phi*nodes_cordinates[:,2];

    # derivative for computation of noradmals
    dy = dphi*nodes_cordinates[:,1];
    dz = dphi*nodes_cordinates[:,2];

    # noradmals in x and y
    ny = -dz;
    nz = dy;

    return y,z,ny,nz;
end

function modifiedshapefunction(xi)

    ## returns the shapefunction and its derivate for linear element
    # xi is a vector of size [n,1]
    # mphi is matrix of size [n,2]

    mphi = [-0.5*ones(size(xi)) 0.5*ones(size(xi))];
    return mphi;
end

## reading calulation file
function readcalcfile(name)
    # Function able to read the input file *.dat, which contains the data of the pre-process.
    # INPUT DATA
    # name --> path and name of the calculation file.


    Nelem = 0; Nnodes=0; Nnodes_IntBEM=0;
    N = readlines(name);

    for i=1:length(N)

        if N[i]=="Number of Elements:"
            Nelem =  parse(Int64,N[i+1])
        end

        if N[i]=="Number of Nodes:"
            Nnodes =  parse(Int64,N[i+1])
        end


        if N[i]=="N.nodes Interior_BE"
            Nnodes_IntBEM =  parse(Int64,N[i+1])
        end

    end
    posnodes = zeros(Nnodes,4);
    conect_mat= zeros(Nnodes,5);
    mat_IntBEM = zeros(Nnodes_IntBEM,3);

    for i=1:length(N)
        if N[i]=="Node X Y Z"
            for j = 1:Nnodes
                posnodes[j,:]=transpose([parse(Float64, ss) for ss in split(N[i+j])])
            end
        end

        if N[i]=="Element    Node(1)   Node(2) Node(3)   Material"
            for j = 1:Nnodes
                conect_mat[j,:]=transpose([parse(Int64, ss) for ss in split(N[i+j])])
            end
        end
        if Nnodes_IntBEM!=0
            if N[i]=="N.nodes Interior_BE"
                for j = 1:Nnodes_IntBEM
                    mat_IntBEM[j,:]=transpose([parse(Int64, ss) for ss in split(N[i+1+j])])
                end
            end
        end
    end

    BGNodes = posnodes;
    BGNodes[:,3] = -BGNodes[:,3];
    posnodes=posnodes[:,2:4];
    conect_mat=conect_mat[:,2:5];

    return posnodes,conect_mat,Nelem,Nnodes,mat_IntBEM,BGNodes;

end

function Connectivity(mat_IntBEM,BGNodes,conect_mat)


    while (isempty(mat_IntBEM)!=1)

        IntCon[1,:] = mat_IntBEM[1,:]; # initialisation

        O = BGNodes[setdiff((conect_mat[IntCon[:,1],1:end-1]),IntCon[:,2:end]),2:end];
        A = BGNodes[(IntCon[:,2]),2:end];
        B = BGNodes[(IntCon[:,3]),2:end];

        c = (A[:,1]-O[:,1])*(B[:,2]-O[:,2])-(A[:,2]-O[:,2])*(B[:,1]-O[:,1]);

        # for the case of closed perimeter in interior case we want c to be
        # greater than zero i.e. counterclockwise

        if c > 0
            mat_IntBEM[:,[2,3]] = mat_IntBEM[:,[3,2]];
        end
        for i=2:length(mat_IntBEM) # connectivity matrix for 1st perimeter

            IntCon[i,:] = mat_IntBEM[mat_IntBEM[:,2] == IntCon[i-1,3],:];
            if  mat_IntBEM[mat_IntBEM[:,2] == IntCon[i,3],:] == IntCon[1,:]
                break;
            else
                continue;
            end
        end

        (~,~,rows) = intersectrows(IntCon,mat_IntBEM); # getting common elements to remove conectivity matrix of 1st perimeter from connity of all perimeters


        #removing rows common in matIntBEM
        mat_IntBEM[setdiff(1:end, rows), :]; # new matIntBEM for next closed perimeter


    end

    return IntCon;
end
function dictionary(IntGConM,BGNodes)
    GConM=IntGConM;

    GConM=GConM[:,2:end];

    NG = GConM[:,1];  # Global Nodes vector

    NL = 1:length(NG);   # Local Nodes vector

    mapG2L = Dict(NG[i] => NL[i] for i = 1:length(NG));

    GNodes=zeros(length(NG),4);

    LConM=zeros(length(NG),3);

    for i = 1:length(NG)
        LConM[i,:] =[i, mapG2L[GConM[i,2]],mapG2L[GConM[i,3]]];

        GNodes[i,:] = BGNodes[GConM[i,2],:];

        NGlob[i] = NG[i];
    end

    return LConM,NGlob,LNodes;
end
using LinearAlgebra
function gaussrule(n)
    ## Returns the weights and base points for the Gauss numerical integration foradmula with 'n' points
    # n is a scalar
    # xi is a row vector of dimension [n,1]
    # wf is a row vector of dimension [n,1]


    # Gauss base points and weight factors calculation taken from
    # a user posted function in the Mathworks website:
    # Concepts on page 93 of
    # 'Methods of Numerical Integration' by
    # Philip Davis and Philip Rabinowitz yield
    # the base points and weight factors.
    #
    #          Howard Wilson
    #          Department of Engineering Mechanics
    #          University of Alabama
    #          Box 870278
    #          Tuscaloosa, Alabama 35487-0278
    #          Phone 205 348-1617
    #          Email address: HWILSON @ UA1VM.UA.EDU


    u = collect(1:n-1)./sqrt.(collect(2 *(1:n-1)).^2 .- 1);
    (xi,vc) = LinearAlgebra.eigen(diagn(u,-1)+diagn(u,1));# think on this
    (xi,k) = sortwithindex((xi));
    wf = 2*vc[1,k]'.^2;

    return xi,wf;
end

## Greens function for acoustics
function AcFullSpace2dot5D(C,kx,y,z,w,ny,nz)
    ## function to comppute greens pressure and displacement matrices

    ## Some possible errors in the input data
    if iscolumn(y)==0;
        error("myApp:argChk","y must be a column vector");
    end
    if iscolumn(z)==0
        error("myApp:argChk","z must be a column vector");
    end
    if iscolumn(kx)==0
        error("myApp:argChk","kx must be a column vector")
    end

    if !(@isdefined ny) && !(@isdefined nz)
        error("myApp:argChk","if ny is an input, nz must be also");
    end

    Nkx = length(kx);  Ny=length(y); Nz = length(z);

    if Ny!= Nz
        error("y and z should have same lengths")
    end

    if isscalar(w)==0 && isscalar(kx)==0
        if length(w)!= length(kx)
            error("w and kx should have same lengths")
        end
    end

    if @isdefined ny
        nyz=[ny nz];
        ny=ny./sqrt.(sum(nyz.^2,dims=2));
        nz=nz./sqrt.(sum(nyz.^2,dims=2));
    end

    # Displacements Green's functions
    # Geometric parameters
    if isempty(C)
        C = 340*(1+2*0.0*1im); # acoustic wave speed
    end

    r = sqrt.(y.^2+z.^2);
    Nr =length(r);

    gy = y./r;
    gz = z./r;

    G = zeros(ComplexF64,1,Nkx,Nr);
    H = zeros(ComplexF64,1,Nkx,Nr);

    for i = 1:length(kx)
        if isscalar(w)==1
            if w == 0.0
                w = 0.001*2*pi;
            end
            ka = sqrt(Complex(kx[i]^2-(w/C)^2));

        else
            ka = sqrt(Complex(kx[i]^2-(w[i]/C)^2));
        end
        Gtemp = (1/(2*pi)).*besselk.(0,ka*r);
        Htemp = (-ka/(2*pi)).*besselk.(1,ka*r).*(gy.*ny+gz.*nz);
        G[:,i,:] = Gtemp;
        H[:,i,:] = Htemp;

    end
    G = reshape(G,(1,1,Nkx*Nr));
    H = reshape(H,(1,1,Nkx*Nr));

    return G,H;
end

function AcFullSpace2dot5DBEM(C,kx,y,z,w)
    ## function to comppute greens pressure and displacement matrices

    ## Some possible errors in the input data
    if iscolumn(y)==0;
        error("myApp:argChk","y must be a column vector");
    end
    if iscolumn(z)==0
        error("myApp:argChk","z must be a column vector");
    end

    Nkx = length(kx);  Ny=length(y); Nz = length(z);

    if Ny!= Nz
        error("y and z should have same lengths")
    end

    if isscalar(w)==0 && isscalar(kx)==0
        if length(w)!= length(kx)
            error("w and kx should have same lengths")
        end
    end


    r = sqrt.(y.^2+z.^2);
    rs = unique(r); Nr =length(rs);

    # Displacements Green's functions
    # Geometric parameters
    if isempty(C)
        C = 340*(1+2*0.0*1im); # acoustic wave speed
    end
    # computing the greesn functions

    G = zeros(ComplexF64,1,Nkx,Nr);
    H = zeros(ComplexF64,1,Nkx,Nr);


    for i = 1:length(kx)
        if isscalar(w)==1
            if w ==0
                w = 0.001*2*pi;
            end
            ka = sqrt(Complex(kx[i]^2-(w/C)^2));
        else
            ka = sqrt(Complex(kx[i]^2-(w[i]/C)^2));
        end
        Gtemp = (1/(2*pi)).*besselk.(0,ka*rs);
        Htemp = (-ka/(2*pi)).*besselk.(1,ka*rs);
        G[:,i,:] = Gtemp;
        H[:,i,:] = Htemp;

    end

    G = reshape(G,(1,1,Nkx*Nr)); #G = G[:,:,1];
    H = reshape(G,(1,1,Nkx*Nr)); #H = H[:,:,1];

    return G,H,rs;
end

function AcFullSpace2dot5DBEMvelocity(H,y,z,r,ny,nz)
    ##function to compute Green"s Velocity

    if @isdefined ny
        nyz=[ny nz];
        ny=ny./sqrt.(sum(nyz.^2,dims=2));
        nz=nz./sqrt.(sum(nyz.^2,dims=2));
    end

    gy = y./r;
    gz = z./r;

    Nr = length(r);

    ang = (gy.*ny+gz.*nz);
    ang = reshape(ang,(1,1,1,Nr));


    H = H.*ang;
    (s1,s2,s3,s4,s5) = size(H)
    Hnew = reshape(H,(1,1,s3*s4));

    return Hnew;
end

## Auxilary functions and macros

function iscolumn(v)
    # function to check if the vector is column
    m=ndims(v)
    if m == 1
        return 1
    elseif m==0
        return 1

    else
        return 0
    end

end

function isscalar(x)
    # function to check if the input is scalar
    if ndims(x) ==1 && length(x) == 1
        return 1
    elseif ndims(x)==0
        return 1
    else
        return 0
    end

end

function ismembertol(a,b,val)
    ## functional equivalent of matlab ismember
    index = zeros(Int,length(a),1)
    for i = 1:length(a)
        index[i]=findfirst(a[i] .== b)

    end

    return val,index;
end

function ismember(a,b)
    ## functional equivalent of matlab ismember
    Lia = in.(a, [b])
    LocB = findall(Lia)

    return Lia, LocB;
end

import Base.@isdefined
macro isdefined(var)
    quote
        try
            local _ = $(esc(var))
            true
        catch err
            isa(err, UndefVarError) ? false : rethrow(err)
        end
    end
end

function intersectrows(ms::Array...)
    t = map(x->Dict(x[2][i,:]=>(x[1],i) for i=1:size(x[2],1)),enumerate(ms))
    u = intersect(map(keys,t)...)
    return (u,map(x->[x[r][2] for r in u],t)...)
end

function mtimesx(A,B,moperandi)

    ## function to compute the matrix multiplication for multidimensional matrices
    ## A and B are two matrices have same dimesnion
    ## moperandi is the variable that defines the how the operation will be perforadmed
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
    elseif ndims(B)-ndims(A)==1
        if ndims(A) == 3 && ndims(B)==4
            A = repeat(A,inner=(1,1,1,1))
        else
            error
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

        if moperandi == "" ## using simd
            C = zeros(ComplexF64,sa1,sb2,sa3);

            for i = 1:sb3
                C[:,:,i] = A[:,:,i]*B[:,:,i];
            end

        elseif moperandi == "threaded" ## using multithreading
            C = zeros(sa1,sb2,sa3);

            Threads.@threads for i = 1:sb3
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

function sortwithindex(n)
    m = sort(n);
    Lia = in.(m, [n])
    index = findall(Lia)
    return m,index;
end

function diagn(m,iden)
    ##  places the elements of vector v on the kth diagonal except main.
    # iden > 0 is above the main diagonal, and k < 0 is below the main diagonal.

    k = abs(iden);
    n = zeros(length(m)+k,length(m)+k)
    if iden > 0
        #for upper diagonal
        for i =1:length(m)+k
            for j = 1:length(m)+k
                if  i<j && j==i+(k)
                    n[i,j] = m[i]
                end
            end
        end
    end
    if iden < 0
        #for lower diagonal
        for i = 1:length(m)+k
            for j = 1:length(m)+k
                if  i > j && i == j+(k)
                    n[i,j] = m[j]
                end
            end
        end
    end
    return n;
end

function findvals(mat,qpoint)

    nd = ndims(mat);
    (ids) = findall(mat .== qpoint);

    nodeid = zeros(Int8, 2, 1);
    elemid = zeros(Int8, 2, 1);

    if nd ==3

        for i = 1:2
            indx = CartesianIndex(ids[i]);
            nodeid[i] = indx[1];
            elemid[i] = indx[3];
        end
    end
    if nd == 2
        for i = 1:2
            indx = CartesianIndex(ids[i]);
            nodeid[i] = indx[1];
            elemid[i] = indx[2];
        end
    end

    return nodeid,elemid;
end
