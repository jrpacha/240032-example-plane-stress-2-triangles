function Fe=internalConstantForcesQuad(nodes,elem,e,th,f)
    % compute internal force terms when the force f is constant
    % for self-weight f=(0,-\ro*g) being \ro the density ang g gravity
    %
    N=2; %number of GaussPoints
    [w,ptGaussRef]=gaussValues2DQuad(N);
    %
    % First compute Ke, Fe for each element
    %    
    v1=nodes(elem(e,1),:);
    v2=nodes(elem(e,2),:);
    v3=nodes(elem(e,3),:);
    v4=nodes(elem(e,4),:); 
    vertices=[v1;v2;v3;v4];
    % Shape functions
    Psi1=@(x,y)(1-x).*(1-y)/4;
    Psi2=@(x,y)(1+x).*(1-y)/4;
    Psi3=@(x,y)(1+x).*(1+y)/4;
    Psi4=@(x,y)(1-x).*(1+y)/4;
    shapeFunctions = @(x,y)[Psi1(x,y),Psi2(x,y),Psi3(x,y),Psi4(x,y)];        
    % Shape function derivatives
    dPsi11=@(x,y) -(1-y)/4;
    dPsi21=@(x,y) (1-y)/4;
    dPsi31=@(x,y) (1+y)/4;
    dPsi41=@(x,y) -(1+y)/4;
    dPsi12=@(x,y) -(1-x)/4;
    dPsi22=@(x,y) -(1+x)/4;
    dPsi32=@(x,y) (1+x)/4;
    dPsi42=@(x,y) (1-x)/4;
    % Derivative matrix 2x4
    Jacob =@(x,y) [dPsi11(x,y), dPsi21(x,y),dPsi31(x,y),dPsi41(x,y);...
                   dPsi12(x,y), dPsi22(x,y),dPsi32(x,y),dPsi42(x,y)];

    % Compute the corresponding Gaussian points on the domain
    % evaluate Shape functions on Gaussian reference points
    xx = ptGaussRef(:,1);
    yy = ptGaussRef(:,2);
    % evaluate Jacobian contribution for each point 
    % We use a Matlab cell variable in order to load each matrix 
    numPtG=size(xx,1);  
    for i=1:numPtG
        Jtilde{i}=inv(Jacob(xx(i),yy(i))*vertices);
        evalDetJacob(i) = det(Jacob(xx(i),yy(i))*vertices);
        Jacobia{i}=Jacob(xx(i),yy(i)); %derivatives of the shape functions
        shapeFun{i}=shapeFunctions(xx(i),yy(i));
    end
    %
    %
    for j=1:4   
        suma=0;
        for ptG=1:numPtG
            suma=suma+w(ptG)*shapeFun{1,ptG}(j)*evalDetJacob(ptG);
        end
        shInt(j)=suma;
    end
    intC=[shInt(1),shInt(1),shInt(2),shInt(2),shInt(3),shInt(3),shInt(4),shInt(4)];

    fx=f(1);
    fy=f(2);
    Fe=th*intC.*[fx,fy,fx,fy,fx,fy,fx,fy];
    Fe=Fe';
