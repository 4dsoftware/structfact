F=imread('C:\Users\Daniel\Reseach\Proj #6 Phase Digram\Data\solid 600 600.tif');
Fd = double(F(:,:,1));
FF = Fd;
b  = bpass(FF,0.3,40);             % 40 for 2.8 um, 8 for 1 um
bb = mat2gray(b);
pk = pkfnd(bb,0.5,20);           % 20 for 2.8 um, 7 for 1 um
cnt = cntrd(bb,pk,27,0);         % 27 for 2.8 um, 9 for 1 um

%Calculate structure factor
sz=size(cnt);
N=sz(1);
qx=-10:0.5:10;
qy=-10:0.5:10;
for k=1:length(qx)
    for l=1:length(qy)
        sum=0;
        for m=1:N
            for n=1:N
                sum=sum+exp(i*((cnt(m,1)-cnt(n,1))*qx(k)+(cnt(m,2)-cnt(n,2))*qy(l)));
            end
        end
        S(k,l)=sum/N;
    end
end

