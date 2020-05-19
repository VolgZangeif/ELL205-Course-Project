function S = SecretImageRetrieval(H,K,O)
%H is the suspected image
%m,n are the dimensions of secret image
%K is the secret key
%O is the ownership share
[s,t] = size(O);
m = s/2; n = t/2;
M = MasterShareCreation(H,ones(m,n),K);
S = GetSecretImage(M,O);
end

function S = GetSecretImage(M,O)
[c,d] = size(M);
m = c/2; n = d/2;
S_ = zeros(c,d);
for k = 1:c
    for l = 1:d
        if(M(k,l) == 255 && O(k,l) == 255) %White-White overlap
           S_(k,l) = 255;
        else
            S_(k,l) = 0;                    %Else result of overlap will be black
        end
    end
end
A = ones(1,m) * 2;
B = ones(1,n) * 2;
C = mat2cell(S_,A,B);
S = zeros(m,n);
for k = 1:m
    for l = 1:n
        if(sum(C{k,l},'all') < 510)
            S(k,l) = 0;
        else
            S(k,l) = 255;
        end
    end
end
end



function [M,O] = MasterShareCreation(hostImage,embedImage,K)
% K is the secret key
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[M,N] = size(hostImage);
[m,n] = size(embedImage);

hostImage = double(hostImage);                  %For SVD
embedImage = imbinarize(embedImage,0.7)*255;    %For making secret image BnW

rng(K);

x = floor(M/4);
y = floor(N/4);

randrows = randi(x,m,n);
randcols = randi(y,m,n);

A = ones(1,x) * 4;
B = ones(1,y) * 4;
hostImage = hostImage(1:(x*4),1:(y*4));
X = mat2cell(hostImage,A,B);   %This cells contains all 4x4 blocks of the image, in a 2D cell
for k = 1:m
    for l = 1:n
        a(:,:,k,l) = X{randrows(k,l),randcols(k,l)};
        a(:,:,k,l) = fft2(a(:,:,k,l));
    end
end

% for i = 1:m
%     for j = 1:n
%         a(:,:,i,j) = fft2(a(:,:,i,j))
%     end
% end

b = SVDStep4(a);

b = BinaryMap(b);

M = MasterShare(b);

O = OwnershipShare(M,embedImage);

% for k = 1:M/4
%     for l = 1:N/4
%         a(:,:,(k-1)*(N/4) + l) = fft2(a(:,:,(k-1)*(N/4) + l))   %linearising the cell into an array
%     end
% end
% a(:,:,:) = fft2(a(:,:,:))

% for i = 1:m
%     for j = 1:n
%         b(i,j) = 
      %size(ans) = 4x4x(M*N/16), access blocks by ans(:,:,i)  
end

function sig = SVDStep4(X)
%Step 4
%X is a 4*4*m*n matrix
[~,~,m,n] = size(X);
sig = ones(m,n);
for i = 1:m
    for j = 1:n
        temp = svd(X(:,:,i,j));
        sig(i,j) = temp(1);
    end
end
end

function B = BinaryMap(X)
%Step 5
av = mean(X,'all');
[m,n] = size(X);
B = zeros(m,n);
for k = 1:m
    for l = 1:n
        if(X(k,l) >= av)
            B(k,l) = 255;
        end
    end
end
end

function M = MasterShare(B)
[m,n] = size(B);
C = cell(m,n);
x = [255 0; 0 255];
y = [0 255; 255 0];
for k = 1:m
    for l = 1:n
        if(B(k,l) == 255)
            C{k,l} = x;
        else
            C{k,l} = y;
        end
    end
end
M = cell2mat(C);
end

function O = OwnershipShare(M,S)
[m,n] = size(S);
A = ones(1,m) * 2;
B = ones(1,n) * 2;
X = mat2cell(M,A,B);  %Converting the master share in a cell of 2x2 blocks
C = cell(m,n);
x = [255 0; 0 255];
y = [0 255; 255 0];
for k = 1:m
    for l = 1:n
        if(S(k,l) == 255 && isequal(X{k,l},x)) 
            C{k,l} = x;
        elseif(S(k,l) == 255 && isequal(X{k,l},y))
            C{k,l} = y;
        elseif(S(k,l) == 0 && isequal(X{k,l},x))
            C{k,l} = y;
        elseif(S(k,l) == 0 && isequal(X{k,l},y))
            C{k,l} = x;
        end
    end
end
O = cell2mat(C);
% O = C;
end
