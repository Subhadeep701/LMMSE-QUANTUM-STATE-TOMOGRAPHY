function [sigma_operator,Sic_POVM] = sigpov(n)
M=eye(n);D=zeros(n,n,n,n);sigma_operator=zeros(n,n,n^2-1);l1=1;
omega=exp(2*pi*1i/n);X=zeros(n,n);
%%%%%% generalized basis operator generation%%%%%%%%%%%
for j1=1:n-1
    for k1=(j1+1):n
        sigma_operator(:,:,l1)=M(:,j1)*M(:,k1)'+M(:,k1)*M(:,j1)';
        l1=l1+1;
        sigma_operator(:,:,l1)=1i*(-M(:,j1)*M(:,k1)'+M(:,k1)*M(:,j1)');
        l1=l1+1;
    end
    temp_op=zeros(n,n);
    for j11=1:j1
    temp_op=M(:,j11)*M(:,j11)'+temp_op;
    end
    sigma_operator(:,:,l1)=sqrt(2/(j1*(j1+1)))*(temp_op-j1*M(:,j1+1)*M(:,j1+1)');
    l1=l1+1;
end
%%%%SIC POVM GENERATION %%%%%%%%%%%%%%%%
for i=1:n
    if(i<n)
      X=M(:,i+1)*M(:,i)'+X;
    else
     X=X+M(:,1)*M(:,i)'; 
    end
end
Z=diag(omega.^(0:n-1));
for i=0:n-1
    for j=0:n-1
        D(:,:,i+1,j+1)=omega^(i*j/2)*(X^j)*(Z^i);
    end
end
fiducial={};
fiducial{1,2}=[3.7637620719571985e-01+1i*2.7760878126176896e+00;
1.1538819157681195e+00+1i*8.7835540934078105e-01];

fiducial{1,3}=[1.9982313502573593e-01-1i*6.3035160323234174e-01;
-2.8913064462318849e-01-1i*3.1990601498724591e-01;
-1.5685828957940073e-01-1i*1.6832169782506001e-01];

fiducial{1,4}=[2.0647049017782240e-01-1i*6.1045242686444368e-01;
-2.1472818276555916e-01+1i*4.8655405452956124e-01;
-2.5285767439082263e-01-1i*8.5522877261474936e-02;
-3.1020136683740951e-01-1i*9.4588221909820591e-01];

fiducial{1,5}=[-1.1142004096427088e-01-1i*6.1880005054111276e-01;
-1.5631170366856961e-02+1i*7.6206148082138669e-01;
3.1404308925728000e-02+1i*2.7378179627692173e-01;
-4.8968143551791687e-01-1i*2.7981552119811817e-02;
-7.3754420009166666e-01-1i*8.6513718163831332e-01];

fiducial{1,6}=[3.2268150227110759e-01-1i*2.5011304940978629e-01;
-1.0079381072262217e-01-1i*7.0341023772541400e-01;
-2.3090135014767987e-01+1i*3.5239719128103686e-01;
-5.0606388793039458e-01-1i*7.8602647353058430e-01;
6.4138470375045400e-02+1i*5.2699200923766354e-01;
6.3737122541549660e-02-1i*6.1877391065498652e-02];

fiducial{1,7}=[2.1335483082419693e-01-1i*1.8302878074462128e-01;
-1.2715199183962100e-01-1i*6.5679038320756633e-01;
-4.2725214851205023e-01-1i*7.7771494980375400e-01;
-3.7253979323384333e-02+1i*2.9153304859443246e-01;
1.8401836082089468e-01+1i*3.9127381792749028e-01;
7.9548090723315257e-02-1i*1.7419950566352835e-01;
-3.6951179931098854e-01+1i*5.4444292952246842e-01];

fiducial{1,8}=[1.9295371721244883e-01-1i*1.5838827120127838e-03;
1.1790390770807807e-01-1i*1.3012384240892941e-01;
-2.9775389568059210e-02-1i*3.2859821728350885e-01;
-6.1507901774529909e-02-1i*3.3325514776353940e-02;
3.1906203109556192e-02+1i*1.0051076498273709e-01;
-4.9748566120997102e-02+1i*7.6305553008150978e-02;
-1.4225233958574424e-01-1i*1.1370278046829319e-01;
-2.0505880353199424e-01-1i*9.5481956657053751e-02];

fiducial{1,9}=[-9.4918510710444076e-02-1i*6.2052177096661457e-02
-4.1763584561678180e-01+1i*4.0112054842789367e-02
1.4916460428113936e-02-1i*2.8166025202824053e-02
2.5836495997042158e-01+1i*1.1934383099886423e-01
-2.5109322409703827e-01+1i*3.0591070729164466e-01
4.8978471703472004e-01+1i*7.2411078337213178e-02
-1.6344644925997509e-01-1i*5.7291653902205406e-02
2.6646584466803946e-01-1i*1.2600778264097862e-01
-1.2127647152891283e-01+1i*7.5423757711034378e-02];

a=fiducial{1,n};
a=a/norm(a);Sic=zeros(n,n^2);k=1;Sic_POVM=zeros(n,n,n^2);
for i=1:n
    for j=1:n
        Sic(:,k)=D(:,:,i,j)*a;
        k=k+1;
    end
end
for l=1:k-1
Sic_POVM(:,:,l)=Sic(:,l)*Sic(:,l)'/n;
end
end

