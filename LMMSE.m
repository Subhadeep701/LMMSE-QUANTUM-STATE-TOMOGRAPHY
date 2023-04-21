clc;
clear all;
close all;
Lin=8;
B={};n=3;sample=16000;clnumber=80;
N=sample/clnumber;
regular_param = 0.0001;
stat=rand(1,Lin);stat=stat/sum(stat);
%stat=[1/3 1/3 1/3];
%stat=[0.9 0.1];
LL=zeros(1,n^2);
%v=zeros(4,n);
%v(1,:)=[ones(1,1),zeros(1,n-1)];v(1,:)=v(1,:)/sum(v(1,:));
%     v(1,:)=[0.45 0.55];
%     v(2,:)=ones(1,n);v(2,:)=v(2,:)/sum(v(2,:));
%     v(3,:)=rand(1,n); v(3,:)=v(3,:)/sum(v(3,:));
%     v(4,:)=[0.45 0.55];v(4,:)=v(4,:)/sum(v(4,:));
%%%%%%%%%%%%%%%%%%%%%%target state%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i2=1:Lin
    X = (randn(n) + 1i*randn(n))/sqrt(2);
    [Q,R] = qr(X);
    v=[ones(1,1),zeros(1,n-1)];v=v/sum(v);
    %v=ones(1,n);v=v/sum(v);
    % v=rand(1,n); v=v/sum(v);
    %v=[ones(1,3),zeros(1,n-3)];v=v/sum(v);
    D1{1,i2}=Q'*(diag(v))*Q
    %D1{1,i2}=Q'*(diag(v(i2,:)))*Q;
end
noise_cov_mtx=zeros(n^2,n^2);
PT_LS=zeros(1,n^2);
%%%%%%%%%%%%%%%%%%%%%% Generalized Pauli and SIC-POVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sigma_operator,tetrahedron_2qb]=sigpov(n);

%%%%%%%%%%%%%%%%H matrix for LINEAR INVERSION%%%%%%%%%%%%%%%%%%%%%%
for d=1:n^2
    for c=1:(n^2)-1
        Anew(d,c)=sum(diag(tetrahedron_2qb(:,:,d)*sigma_operator(:,:,c)))*sqrt(n/(2*(n-1)));
    end
end
%U=[0.01 0.1 0.5];
U=0;
for r=1:length(U)
    E_thet_cov=zeros(n^2-1,n^2-1);noise_ls_cov=zeros(n^2,n^2);
    E_thet=zeros(n^2-1,1);E_Y=zeros(n^2,1);E_thet_noise1=zeros(n^2-1,n^2);E_noise_thet1=zeros(n^2,n^2-1);
    LL=zeros(1,n^2);R_xy=zeros(n^2-1,n^2);R_yy=zeros(n^2,n^2);
    U1=U(r);
    for i=1:clnumber
        PT=zeros(1,n^2);
        C=zeros(1,n^2);
        l=measure(stat);
        D=D1{1,l};
        STORE_D{1,i}=D;
        WD=(1-U1)*D+U1*eye(n,n)/n;
        probs_2qbe = qmt(WD, tetrahedron_2qb);
        counts_2qbe = histc(rand(N,1), [0; cumsum(probs_2qbe)]);
        counts_2qbe = counts_2qbe(1:end-1);
        PT=counts_2qbe/N;

        PTN{1,i}=PT;
        Y(:,i)= (PT-1/n)*(n/(n-1));
        thet{1,i}=inv(Anew'*Anew)*Anew'*Y(:,i);
        E_Y=Y(:,i)+E_Y;
        E_thet=thet{1,i}+E_thet;
    end
    E_Y=E_Y/clnumber;
    E_thet=E_thet/clnumber;
    for i11=1:clnumber
        E_thet_cov=(thet{1,i11}-E_thet)*(thet{1,i11}-E_thet)'+E_thet_cov;
    end
    E_thet_cov=E_thet_cov/clnumber;
    noise_ls_vec=Y-E_Y;
    for n1=1:clnumber
        noise_ls_cov1=noise_ls_vec(:,n1)*noise_ls_vec(:,n1)';
        noise_ls_cov=noise_ls_cov1+noise_ls_cov;
        E_thet_noise1=(thet{n1}-E_thet)*noise_ls_vec(:,n1)'+ E_thet_noise1;
        E_noise_thet1=noise_ls_vec(:,n1)*(thet{n1}-E_thet)'+E_noise_thet1;
        R_xy=(thet{n1}-E_thet)*(Y(:,n1)-E_Y)'+R_xy;
        R_yy=(Y(:,n1)-E_Y)*(Y(:,n1)-E_Y)'+R_yy;
    end
    noise_ls_cov=noise_ls_cov/clnumber;
    E_thet_noise1=E_thet_noise1/clnumber;
    E_noise_thet1=E_noise_thet1/clnumber;
    R_xy=R_xy/clnumber;
    R_yy=R_yy/clnumber;
    R_yy=R_yy+0.000001*eye(n^2)
    %     R_xy=E_thet_cov*Anew'+E_thet_noise1;
    % R_yy=(Anew*E_thet_cov*Anew'+Anew*E_thet_noise1+E_noise_thet1*Anew'+noise_ls_cov);

    %%%%%Testing%%%%%
    %     for j2=1:Lin
    %      for i=1:n^2-1
    %      thet_check(i,j2)=sum(diag(D1{1,j2}*sigma_operator(:,:,i)))*sqrt(n/(2*(n-1)));
    %      end
    %     end
    %     E_thet_check=0.1*thet_check(:,1)+0.3*thet_check(:,2)+0.6*thet_check(:,3);
    %     thet_check_diff=thet_check-E_thet_check;
    %     E_thetcheck_cov=zeros(n^2-1,n^2-1);
    %     for j3=1:Lin
    %         E_thetcheck_cov=stat(j3)*thet_check_diff(:,j3)*thet_check_diff(:,j3)'+E_thetcheck_cov;
    %     end
    %     Dln_LS1= eye(n,n)/n;
    %     for e=1:n^2-1
    %         Dln_LS1=thet_check(e,3)*sigma_operator(:,:,e)*sqrt((n-1)/(2*n))+ Dln_LS1;
    %     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%NEW QUANTUM STATE%%%%%%
    repe =20;
    for k1=1:repe
        Cnew=zeros(1,n^2);
        samplenew=200;
        llast=measure(stat);
        Dlast=D1{1,llast};

        WDN=(1-U1)*Dlast+U1*eye(n,n)/n;

        probs_2qb = qmt(WDN, tetrahedron_2qb);


        counts_2qb = histc(rand(samplenew,1), [0; cumsum(probs_2qb)]);
        counts_2qb = counts_2qb(1:end-1);
        PTnew=counts_2qb/samplenew;
        Y_new=(PTnew-1/n)*(n/(n-1));
        thetnew=pinv(Anew)*Y_new;
        noisecov_2=diag((PTnew-PTnew.^2));
        %Lmmse_thet =E_thet+ inv(inv(E_thet_cov) + Anew'*inv(noise_ls_cov)*Anew)*Anew'*inv(noise_ls_cov)*(Y_new-Anew*E_thet);
        %Lmmse_thet=(E_thet)+ E_thet_cov*Anew'/(Anew*E_thet_cov*Anew'+ noise_ls_cov +noisecov_2)*(Y_new-Anew*E_thet);
        Lmmse_thet=R_xy*inv(R_yy)*(Y_new-E_Y)+E_thet;
        % Lmmse_thet=E_thet+E_thet_cov*Anew'*inv(Anew*E_thet_cov*Anew'+noise_ls_cov + regular_param*eye(n^2))*(Y_new-Anew*E_thet);
        Dln= eye(n,n)/n;
        Dln_LS= eye(n,n)/n;
        for e=1:n^2-1
            Dln=Lmmse_thet(e)*sigma_operator(:,:,e)*sqrt((n-1)/(2*n))+Dln;
            Dln_LS=thetnew(e)*sigma_operator(:,:,e)*sqrt((n-1)/(2*n))+ Dln_LS;
        end
        Dln=proj_spectrahedron(Dln);
        Dln_LS=proj_spectrahedron(Dln_LS);
        T_ls(r,k1)=sum(abs(eig(Dlast-Dln_LS)))/2;
        T_LMMSE(r,k1)=sum(abs(eig(Dlast-Dln)))/2;
    end
end
co_ls=sum(T_ls,2)/repe;
co_lmmse=sum(T_LMMSE,2)/repe;
for j2=1:Lin
     for i=1:n^2-1
         thet_check(i,j2)=sum(diag(D1{1,j2}*sigma_operator(:,:,i)))*sqrt(n/(2*(n-1)));
     end
     for i1=1:n^2
         PT_check(i1,j2)=sum(diag(D1{1,j2}*tetrahedron_2qb(:,:,i1)));
     end
 end
 Y_check=(PT_check-1/n)*(n/(n-1));
 E_thet_check=zeros(n^2-1,1);
 E_Y_check=zeros(n^2,1);
 E_thetcheck_cov=zeros(n^2-1,n^2-1);
 Y_cov_check=zeros(n^2,n^2);
 R_xy_check=zeros(n^2-1,n^2);
 for j3=1:Lin
     E_thet_check=stat(j3)*thet_check(:,j3)+E_thet_check;
     E_Y_check=stat(j3)*Y_check(:,j3)+E_Y_check;
     thet_check_diff=thet_check-E_thet_check;
     Y_check_diff=Y_check-E_Y_check;
  end
 for j3=1:Lin
     E_thetcheck_cov=stat(j3)*thet_check_diff(:,j3)*thet_check_diff(:,j3)'+E_thetcheck_cov;
     Y_cov_check=stat(j3)*Y_check_diff(:,j3)*Y_check_diff(:,j3)'+Y_cov_check;
     R_xy_check=stat(j3)*thet_check_diff(:,j3)*Y_check_diff(:,j3)'+R_xy_check;
 end
 norm(R_yy-Y_cov_check)