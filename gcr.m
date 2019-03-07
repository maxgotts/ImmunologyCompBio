function [pvec,maxaffinvec,affinvec,globjumpvec] = gcr(AgNum)
%This function simulates the GCR.  Typical application:
%
%[pvec,maxaffinvec,affinvec,globjumpvec] = gcr(15);
%
load antigen;
Ag = Ag(AgNum,:);
%Npop = 5e4;
Npop = 2e3;
naiveAb = unidrnd(4,Npop,30);
affin = zeros(Npop,1);
for i = 1:Npop,
    affin(i) = exp(-sum(abs(Ag-naiveAb(i,:)))/30);
end
newAb = naiveAb;
pvec = [];
maxaffinvec = [];
affinvec = [];
maxaffinvec(1) = max(affin);
globjumpvec = [];

[a,b] = hist(affin,100);
figure(1);
clf;
bar(b,a/sum(a));
grid;
xlabel('Affinity');
ylabel('Fraction of the Ab population');
title(strcat('Naive Generation Maximum Affinity  = ',32,32,num2str(max(affin))));

threshold = 0.95;
maxtime = 30;
t = 1;

clc;
while t<=maxtime & max(affin)<threshold,
    p = input('What value of p do you want to use?  ');
    if isempty(p),
        disp('Incorrect key stroke');
        p = input('What value of p do you want to use?  ');
    else end
    if p>1 | p<0,
        disp('p must be between 0 and 1');
        p = input('What value of p do you want to use?  ');
    else end
    tic;
    pvec = [pvec;p];
    H2 = binornd(1,p,Npop,1);
    globjumpind = 0;  
    for i = 1:Npop,
        if H2(i)==1,
            candaffin = [];
            cand = [];
            for j = 1:30,                
                for k = 1:3,
                    candAb = newAb(i,:);
                    candAb(j) = rem(candAb(j)+k,4)+1;
                    cand = [cand;candAb];
                    candaffin = [candaffin;exp(-sum(abs(Ag-candAb))/30)];
                end
            end
            [dummy,Imaxcandaffin] = max(candaffin);
            newAb(i,:) = cand(Imaxcandaffin,:);
        else
            H4 = unidrnd(3,1,30)-2;
            newAb(i,:) = rem(newAb(i,:)+H4,4)+1;
            globjumpind = globjumpind + 1;  
        end
        affin(i) = exp(-sum(abs(Ag-newAb(i,:)))/30);
    end
    globjumpvec = [globjumpvec;globjumpind];
    t = t+1;
    affinvec = [affinvec,affin];
    maxaffinvec = [maxaffinvec;max(affin)];
    disp([pvec maxaffinvec(2:end)]);
    disp('');
    [a,b] = hist(affin,100);
    toc;
    figure(1);
    bar(b,a/sum(a));
    grid;
    xlabel('Affinity');
    ylabel('Fraction of the Ab population');
    title(strcat('Generation ',32,num2str(t),32,'Maximum Affinity  = ',32,32,num2str(max(affin))));
end
figure(2);
clf;
subplot(2,1,1),
plot(max(affinvec),'b.');
hold on
plot(mean(affinvec),'r.');
plot(min(affinvec),'k.');
plot(median(affinvec),'g.');
hold off;
grid;
subplot(2,1,2),
bar(globjumpvec);
Rmean = corrcoef(globjumpvec,mean(affinvec));
Rmedian = corrcoef(globjumpvec,median(affinvec));
Rmax = corrcoef(globjumpvec,max(affinvec));
Rmin = corrcoef(globjumpvec,min(affinvec));
disp(strcat('Correlation between Number of Global Jumps and Mean Affinity =',32,num2str(Rmean(1,2))));
disp(strcat('Correlation between Number of Global Jumps and Median Affinity =',32,num2str(Rmedian(1,2))));
disp(strcat('Correlation between Number of Global Jumps and Max Affinity =',32,num2str(Rmax(1,2))));
disp(strcat('Correlation between Number of Global Jumps and Min Affinity =',32,num2str(Rmin(1,2))));


    