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
# Nworkers needs to be integer specifying the number of workers on parallel computing toolbox
# parallel needs to be 1 for computing with parallel tool box.
# c is wavespeed c == [] implies c = 340 m/s without damping
## The solver function

function AcBEMsolver(route,c,n,w,kx,Ub,FPnt,boundarypressure,BEdomain,sampling)
    ## The boundary Element solver


    if isempty(route)==1 # testting purpose

        ConM = Opts.ConM;  Nodes = Opts.Nodes;

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

    Pf = zeros(Neval,Nkx,Nw);

    if boundarypressure == "yes"

        Pb = zeros(Nnodes,Nkx,Nw);

    else

        Pb = [];

    end


    for nw = 1:Nw
        if nw==2
            tic;
        end

        if sampling=="moving"
            # computing boundary matrices for moving source
            (matHb,matGb,matHf,matGf) = acBEBndMatrices(c,n,w[nw,:],kx[nw,:],xi,wf,phi,Y,Z,Ny,Nz,Vs0,jacobi,nodes,LPnt,FPnt,ConM,Nodes,Vs);

            Vb = permute(repeat(permute(w[nw,:].^2,3,2,1),inner=(size(Ub,1)),1,1).*Ub(nw,:,:),3,2,1);

        else
            # computing boundary matrices for non-moving source
            (matHb,matGb,matHf,matGf) = acBEBndMatrices(c,n,w[nw],kx[nw,:],xi,wf,phi,Y,Z,Ny,Nz,Vs0,jacobi,nodes,LPnt,FPnt,ConM,Nodes,Vs);

            Vb = permute(repeat(permute(w[nw].^2,3,2,1),inner=(size(Ub,1)),1,1).*Ub(nw,:,:),3,2,1);
        end

        rho = 1.24;
        Vb = (1i*1i)*rho*Vb;

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

    # Daigonal Term for acoustics
    (Hs) = AcDiagTerm(c,n,wf,Y,Z,Ny,Nz,jacobi,LPnt,ConM,ElemId,BEdomain);

    # Relative distances for static velocities
    (~,ys,zs,~,nys,nzs) = singDist([],Nodes,ConM,ElemId,LPnt,xi);

    # Static velocities
    (~,Hs0) = AcFullSpace2dot5D(c,0,ys,zs,0,nys,nzs);

    return xi,wf,phi,Y,Z,Ny,Nz,jacobi,nodes,Hs,Hs0,LPnt;
end

function AcDiagTerm(c,n,wf,Y,Z,Ny,Nz,jacobi,LPnt,LConM,ElemId,BEdomain)
    ## function for computing daigonal term for acoustics

    wf = repeat(wf',inner=(1,1,size(LConM,1)));
    He = zeros(size(LPnt,1),1,size(LConM,1));


    (yrel,zrel,~,rrel,ny,nz) = relativedistance([],LPnt,Y,Z,Ny,Nz,jacobi);
    (~,H) = AcFullSpace2dot5D(c,0,yrel,zrel,0,ny,nz); #computing greens function for all relative distance

    for i = 1:size(LPnt,1)

        indx=[size(LConM,1)*n*(i-1)+([1:size(LConM,1)]*n)];

        rm = rrel(indx);

        (~,ibs) = ismembertol(rm,rrel,1e-6);

        Hs4D = reshape(H,1,1,1,length(rrel));
        H0 = reshape(Hs4D[:,:,:,ibs],1,1,size(LConM,1)*n);
        H0 = real(H0);   # real part of static function

        H0r = reshape[H0(1,1,:),n,1,size(LConM,1)] ;
        Hre = zeros(1,1,size(LConM,1)) ;

        Hre[:,:,:] = mtimesx(wf,(H0r[:,:,:].*jacobi[:,:,:])) ;

        for j=1:size(ElemId,1)

            Hre[:,:,ElemId[j,i]] = 0;

        end

        He[i,:,:] = Hre;

    end

    Hs = sum.(He[:,:,:],3);

    if BEdomain !="interior"
        I = repeat(diag(1),inner(size(He,1),1));
        Hs = I-Hs;
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

    # Getting the pressure at the field Points
    Pf = -(MatHf*Pb-MatGf*Vb);

    return  Pf,Pb;

end

function AcIntRepEqn(Pn,Vn,LConM,LNodes,FPnt,nkx)
    ## function to assemble the Green's velocities and pressure matrices for fieldpoint

    MatG = zeros(size(FPnt,1),size(LNodes,1),1,nkx);
    MatH = zeros(size(FPnt,1),size(LNodes,1),1,nkx);

    Pn = reshape(Pn,size(FPnt,1),2,size(LConM,1),nkx);
    Vn = reshape(Vn,size(FPnt,1),2,size(LConM,1),nkx);


    # Assembling the Green's velocities and pressure Matrices
    for i = 1:size(LConM,1)

        MatG[:,[LConM[i,2] LConM[i,3]],:,:] = MatG[:,[LConM[i,2] LConM[i,3]],:,:]+Pn[:,:,i,:];
        MatH[:,[LConM[i,2] LConM[i,3]],:,:] = MatH[:,[LConM[i,2] LConM[i,3]],:,:]+Vn[:,:,i,:];

    end

    MatG=MatG(:,:,:);
    MatH=MatH(:,:,:);

    return MatG,MatH;
end

function AcIntRepEqnMat(c,n,w,kx,wf,phi,Y,Z,Ny,Nz,jacobi,FPnt,LConM)

    ## function to make Green's velocities and Dispalacement  matrices for fieldpoint

    nkx = length(kx); kx= kx'; w=w'; # f and kx must be a column vector
    Jacobi = repeat(jacobi,1,1,1,nkx);

    wf = repeat(wf',inner=(1,1,size(LConM,1)));
    phi = repeat(phi,inner=(1,1,size(LConM,1),nkx));

    Vn = zeros(size(FPnt,1),2,size(LConM,1)*nkx);
    Pn = zeros(size(FPnt,1),2,size(LConM,1)*nkx);

    (yrel,zrel,~,rrel,ny,nz) = relativedistance([],FPnt,Y,Z,Ny,Nz,jacobi);
    (Ps,Us,rs) = AcFullSpace2dot5DBEM(c,kx,yrel,zrel,w);

    for i = 1:size(FPnt,1)

        indx = [size(LConM,1)*n*(i-1)+collect(1:size(LConM,1)*n)];

        rm = rrel(indx);

        (~,ibs) = ismembertol(rm,rs,1e-6);

        Ps4D = reshape(Ps,1,1,nkx,length(rs));
        Vs4D = reshape(Us,1,1,nkx,length(rs));

        P = reshape(Ps4D[:,:,:,ibs],1,1,nkx*size(LConM,1)*n);
        U = AcFullSpace2dot5DBEMvelocity(Vs4D[:,:,:,ibs],yrels,zrels,rm,nys,nzs);


        Pr = reshape(P[1,1,:],nkx,n,1,size(LConM,1));
        Vr = reshape(U[1,1,:],nkx,n,1,size(LConM,1));

        Pr = permute(Pr,[2 3 4 1]);
        Vr = permute(Vr,[2 3 4 1]);

        # vectors for greens pressure on elements
        Pr1 = zeros(1,1,size(LConM,1),nkx);
        Pr2 = zeros(1,1,size(LConM,1),nkx);

        # vectors for greens velocities on nodes of element
        Vr1 = zeros(1,1,size(LConM,1),nkx);
        Vr2 = zeros(1,1,size(LConM,1),nkx);

        # Greens displacemnet matrix element computed at first nodal point for all elements
        Pr1[:,:,:,:] = mtimesx(wf,(phi[:,1,:,:].*Pr[:,:,:,:].*Jacobi[:,:,:,:]));
        # Greens displacemnet matrix element computed at second nodal point for all elements
        Pr2[:,:,:,:] = mtimesx(wf,(phi[:,2,:,:].*Pr[:,:,:,:].*Jacobi[:,:,:,:]));

        # Greens velocities matrix element computed at first nodal point for all elements
        Vr1[:,:,:,:] = mtimesx(wf,(phi[:,1,:,:].*Vr[:,:,:,:].*Jacobi[:,:,:,:]));
        # Greens velocities matrix element computed at second nodal point for all elements
        Vr2[:,:,:,:] = mtimesx(wf,(phi[:,2,:,:].*Vr[:,:,:,:].*Jacobi[:,:,:,:]));

        # converting pressure in 3d matrix for speed
        Pr1=Pr1[:,:,:];
        Pr2=Pr2[:,:,:];

        # converting velocities in 3d matrix for speed
        Vr1=Vr1[:,:,:];
        Vr2=Vr2[:,:,:];

        # the final matrices of pressure and pressure
        Pn[i,:,:] = [Pr1  Pr2];

        Vn[i,:,:] = [Vr1  Vr2];
    end

    return  Pn, Vn;
end

## Boundary Solution Matrices

function acBEBndMatrices(c,n,w,kx,xi,wf,phi,Y,Z,Ny,Nz,Vs0,jacobi,nodes,LPnt,FPnt,LConM,LNodes,Vs)

    ## function to make BE matrices

    # Elements for BE matrices
    (Vn,Pn,CPntId,Lpnt) = AcBndIntEqnMat(c,n,w,kx,xi,wf,phi,Y,Z,Ny,Nz,Vs0,jacobi,nodes,LPnt,LConM,LNodes);

    #Size of kx vector
    nkx = length(kx);

    # Boundary Matrices for Greens velocities and pressure on Boundary
    (MatGb,MatHb) = AcBndIntEqnInt(Vs,Pn,Vn,Lpnt,CPntId,LNodes,LConM,nkx);

    if !isempty(FPnt)
        # Getting Greens velocities and pressure for fieldpoint
        (Pn,Vn)=AcIntRepEqnMat(c,n,w,kx,wf,phi,Y,Z,Ny,Nz,jacobi,FPnt,LConM);

        # Building the velocities and pressure matrices for fieldpoint
        (MatGf,MatHf)=AcIntRepEqn(Pn,Vn,LConM,LNodes,FPnt,nkx);

    else

        MatHf=[]; MatGf=[];

    end

    return MatHb,MatGb,MatHf,MatGf;
end

function AcBndIntEqnInt(Vs,Pn,Vn,Lchf,CPntId,LNodes,LConM,nkx)
    ## function to assemble the Greens velocities and pressure matrices on boundray

    MatG = zeros(Lchf,size(LNodes,1),1,nkx);
    MatH = zeros(Lchf,size(LNodes,1),1,nkx);

    Pn = reshape(Pn,Lchf,2,size(LConM,1),nkx);
    Vn = reshape(Vn,Lchf,2,size(LConM,1),nkx);

    # Assembling the Green's velocities and pressure Matrices
    MatI = zeros(Lchf,size(LNodes,1));

    for i = 1:size(LConM,1)

        MatG[:,[LConM[i,2] LConM[i,3]],:,:] = MatG[:,[LConM[i,2] LConM[i,3]],:,:]+Pn[:,:,i,:];
        MatH[:,[LConM[i,2] LConM[i,3]],:,:] = MatH[:,[LConM[i,2] LConM[i,3]],:,:]+Vn[:,:,i,:];

    end

    for i = 1:size(LNodes,1)

        MatI[i,CPntId[i]] = MatI[i,CPntId[i]]+Vs[i,:];

    end

    MatH=MatH[:,:,:];
    MatG=MatG[:,:,:];
    MatH = MatH+repeat(MatI,inner=(1,1,nkx));

    return MatG,MatH;
end


function AcBndIntEqnMat(c,n,w,kx,xi,wf,phi,Y,Z,Ny,Nz,Us0,jacobi,nodes,LPnt,LConM,LNodes)
    ## function to make the Greens velocities and pressure matrices on boundary

    nkx = length(kx); kx = kx'; w=w'; # f and kx needs to be a column vector
    Jacobi=repeat(jacobi,1,1,1,nkx);

    wf = repeat(wf',inner=(1,1,size(LConM,1)))
    phi = repeat(phi,inner=(1,1,size(LConM,1),nkx));

    (~,ElemId,CPntId) = singLocate(LPnt,LConM,LNodes,nodes);
    l = size(LPnt,1);

    Lpnt = size(LPnt,1);

    Pn = zeros(Lpnt,2,size(LConM,1)*nkx);
    Vn = zeros(Lpnt,2,size(LConM,1)*nkx);

    (yrel,zrel,~,rrel,ny,nz) = relativedistance([],LPnt,Y,Z,Ny,Nz,jacobi);
    (G,H,rs) = AcFullSpace2dot5DBEM(c,kx,yrel,zrel,w);

    for i = 1:size(LPnt,1)

        indx = (size(LConM,1)*n*(i-1)+collect(1:size(LConM,1)*n));

        rm = rrel[indx];

        yrels = yrel[indx];
        zrels = zrel[indx];

        nys = ny[indx];
        nzs = nz[indx];

        (~,ibs) = ismembertol(rm,rs,1e-6);

        Ps4D = reshape(G,1, 1, nkx, length(rs));
        Vs4D = reshape(H,1, 1, nkx, length(rs));

        P = reshape(Ps4D[:,:,:,ibs],1, 1, nkx*size(LConM,1)*n);
        V = AcFullSpace2dot5DBEMvelocity(Vs4D[:,:,:,ibs],yrels,zrels,rm,nys,nzs);

        # reshaping greens pressure vector
        Pr = reshape(P[1,1,:],nkx,n,1,size(LConM,1)) ;
        Pr = permute(Pr,[2 3 4 1]);

        # reshaping greens velocities vector
        Vr = reshape(V[1,1,:],nkx,n,1,size(LConM,1)) ;
        Vr = permute(Vr,[2 3 4 1]);

        # vectors for greens pressure on nodes of element
        Pr1 = zeros(1,1,size(LConM,1),nkx) ;
        Pr2 = zeros(1,1,size(LConM,1),nkx) ;

        # vectors for greens velocities on nodes of element
        Vr1 = zeros(1,1,size(LConM,1),nkx) ;
        Vr2 = zeros(1,1,size(LConM,1),nkx) ;

        ## Regular Integration

        # Greens displacemnet matrix element computed at first nodal point for all elements
        Pr1[:,:,:,:] = mtimesx(wf,(phi[:,1,:,:].*Pr[:,:,:,:].*Jacobi[:,:,:,:]));
        # Greens displacemnet matrix element computed at second nodal point for all elements
        Pr2[:,:,:,:] = mtimesx(wf,(phi[:,2,:,:].*Pr[:,:,:,:].*Jacobi[:,:,:,:]));

        # Greens velocities matrix element computed at first nodal point for all elements
        Vr1[:,:,:,:] = mtimesx(wf,(phi[:,1,:,:].*Vr[:,:,:,:].*Jacobi[:,:,:,:]));
        # Greens velocities matrix element computed at second nodal point for all elements
        Vr2[:,:,:,:] = mtimesx(wf,(phi[:,2,:,:].*Vr[:,:,:,:].*Jacobi[:,:,:,:]));

        ## singular Integration
        if i<l+1

            for j=1:size(ElemId,1)

                (mphi) = modifiedshapefunction(xi);
                mphi = repeat(mphi,1,1,1,nkx);

                k=2*(i-1)+j;

                V0ss = Us0[:,:,collect((k-1)*n+1:k*n)];
                Vss = real(V0ss); # real part for static function
                Vms = real(V0ss); # real part for modified static function

                Vssr = reshape(Vss(1,1,:),[n,1,1]) ;
                Vssr = repeat(Vssr,1,1,1,nkx);

                Vmsr = reshape(Vms[1,1,:],[n,1,1]) ;
                Vmsr = repeat(Vmsr,1,1,1,nkx);

                Vdr1 = mtimesx(wf[:,:,ElemId[j,i],:],(phi[:,1,ElemId[j,i],:].*(Vr[:,:,ElemId[j,i],:]-Vssr).*Jacobi[:,:,ElemId[j,i],:])) ;
                Vdr2 = mtimesx(wf[:,:,ElemId[j,i],:],(phi[:,2,ElemId[j,i],:].*(Vr[:,:,ElemId[j,i],:]-Vssr).*Jacobi[:,:,ElemId[j,i],:])) ;

                V0sr1 = mtimesx(wf[:,:,ElemId[j,i],:],(mphi[:,1,:,:].*(Vmsr).*Jacobi[:,:,ElemId[j,i],:])) ;
                V0sr2 = mtimesx(wf[:,:,ElemId[j,i],:],(mphi[:,2,:,:].*(Vmsr).*Jacobi[:,:,ElemId[j,i],:])) ;

                # singular Greens velocities
                Vr1[:,:,ElemId[j,i],:] = Vdr1+V0sr1 ;
                Vr2[:,:,ElemId[j,i],:] = Vdr2+V0sr2 ;
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

    CPntId = zeros(1,size(LNodes,1));
    NdeId = zeros(2,size(LNodes,1)); # node id of singularity
    ElemId = zeros(2,size(LConM,1)); # element id of singularity

    for i = 1: size(LPnt,1)

        # finding singular element and nodal location of singularity on that element
        Pxy = ones(2,1)*LPnt[i,:];
        pntmat = repeat(Pxy,1,1,size(LConM,1));

        testpnt = sum(abs((nodes[:,:,:]-pntmat)),2);
        (nodeid,elmtid) = findall(testpnt == 0);

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

    nys = zeros(2*length(xi),size(LPnt,1)); # normals on singular element
    nzs = zeros(2*length(xi),size(LPnt,1)); # normals on singular element

    for i= 1:size(LPnt,1)

        y_sing = zeros(2,length(xi));
        z_sing = zeros(2,length(xi));

        ny_sing = zeros(2,length(xi));
        nz_sing = zeros(2,length(xi));

        for j=1:size(ElemId,1)

            NCordinates = LNodes[LConM[ElemId[j,i],2:3],2:end];
            (y,z,ny,nz) = isoparametric(NCordinates,xi);

            y = y-LPnt[i,1];

            if Method!="hlfspc"
                z = z-LPnt[i,2];
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


    for i=1:size(LPnt)

        y[:,i] = Y[:]-LPnt[i,1];

        if Method =="hlfspc"

            z[:,i] = Z[:];

        else

            z[:,i] = Z[:]-LPnt[i,2];
        end
    end

    ny = Ny./jacobi; ny = ny[:];
    nz = Nz./jacobi; nz = nz[:];

    ny = repeat(ny,[size(LPnt,1),1]);
    nz = repeat(nz,[size(LPnt,1),1]);

    YZrel = round.([y[:],z[:]],5);

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
    (X,Z,Nx,Nz,jacobi,LPnt,nodes) = geometrytransform(LConM,LNodes,xi);

    return xi,wf,phi,X,Z,Nx,Nz,jacobi,LPnt,nodes;
end

function geometrytransform(LConM,LNodes,xi)
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

    phi = [0.5*(1-xi) 0.5*(1+xi)];
    dphi = [-0.5*ones(size(xi)) 0.5*ones(size(xi))];
    return phi,dphi;
end

function isoparametric(nodes_cordinates,xi)

    ## Returns cordinates y,z; shape function phi; normals nx,ny; and jacobi as a function of xi
    # nodes_cordinates is matrix of size [2,2] the columns of which represent cordinates x,y
    # xi is vector of size[n,1]
    # y,z,ny,nz,jacobi are vectors of size[n,1]
    # phi is matrix of size[n,2]

    # Linear Shape functions and their derivatives
    (phi,dphi) = shapefunction(xi);

    # Local cordinate of x and y
    y = phi*nodes_cordinates[:,1];
    z = phi*nodes_cordinates[:,2];

    # derivative for computation of normals
    dy = dphi*nodes_cordinates[:,1];
    dz = dphi*nodes_cordinates[:,2];

    # normals in x and y
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

function [xi,wf]=gaussrule(n)

    ## Returns the weights and base points for the Gauss numerical integration formula with 'n' points
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


    u = (1:n-1)./sqrt.((2 *(1:n-1)).^2 .- 1);
    [vc,xi] = eig(diag(u,-1)+diag(u,1));# think on this
    [xi,k] = sort(diag(xi));
    wf = 2*vc[1,k]'.^2;

end

function ismember(a,b)
    ## functional equivalent of matlab ismember
    Lia = in.(a, [b])
    LocB = findall(Lia)

    return Lia, LocB;
end
