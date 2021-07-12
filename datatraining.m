clc;clear;close all;
%membaca gambar
image_folder = 'D:\Project_pcd\data training';
filenames = dir(fullfile(image_folder, '*.jpg'));
total_images = numel(filenames);

for n = 1:total_images
  full_name= fullfile(image_folder, filenames(n).name);
  our_images = imread(full_name);
  img=rgb2gray(our_images);
  img=histeq(img);
%kondisi jika nilai matriks ada 0 semua di tambah 1
if(min(min(img)) == 0)
    img = img + 1;
end
%nilai max pada gambar
maximg = max(max(img));
%panjang dan lebar gambar(pixel)
[panjang,lebar]=size(img);
%membuat matriks kosong sejumlah 4
matrikskosong=zeros(maximg,maximg,4);

%glcm 0 derajat
for k =1:panjang
    for l=1:lebar-1
        matriksbantu0 = img(k,l:l+1);
        matrikskosong(matriksbantu0(1,2),matriksbantu0(1,1),1)=matrikskosong(matriksbantu0(1,2),matriksbantu0(1,1),1)+1;
    end
end
matrikstranspose0 = transpose(matrikskosong(:,:,1));
hasil0 = matrikskosong(:,:,1) + matrikstranspose0;

%glcm 45 derajat
for i =1:panjang-1
    for j=1:lebar-1
        matriksbantu45 = img(i:i+1,j:j+1);
        matrikskosong(matriksbantu45(1,2),matriksbantu45(2,1),2)=matrikskosong(matriksbantu45(1,2),matriksbantu45(2,1),2)+1;
    end
end
matrikstranspose45 = transpose(matrikskosong(:,:,2));
hasil45 = matrikskosong(:,:,2) + matrikstranspose45;

%glcm 90 derajat
for i =1:panjang-1
    for j=1:lebar
        matriksbantu90 = img(i:i+1,j);
        matrikskosong(matriksbantu90(1,1),matriksbantu90(2,1),3)=matrikskosong(matriksbantu90(1,1),matriksbantu90(2,1),3)+1;
    end
end
matrikstranspose90 = transpose(matrikskosong(:,:,3));
hasil90 = matrikskosong(:,:,3) + matrikstranspose90;

%glcm 135 derajat
for i =1:panjang-1
    for j=1:lebar-1
        matriksbantu135 = img(i:i+1,j:j+1);
        matrikskosong(matriksbantu135(1,1),matriksbantu135(2,2),4)=matrikskosong(matriksbantu135(1,1),matriksbantu135(2,2),4)+1;
    end
end
matrikstranspose135 = transpose(matrikskosong(:,:,4));
hasil135 = matrikskosong(:,:,4) + matrikstranspose135;

%normalisasi
normalisasi0=hasil0/sum(sum(hasil0));
normalisasi45=hasil45/sum(sum(hasil45));
normalisasi90=hasil90/sum(sum(hasil90));
normalisasi135=hasil135/sum(sum(hasil135));
%panjang dan lebar gambar yang telah dinormalisasi
[ukuranx,ukurany]=size(normalisasi0);

%kontrast 0 derajat
kontrast0=0;
for i=1:ukuranx
    for j=1:ukurany
        kontrast0 = ((i-j)^2*normalisasi0(i,j))+kontrast0;
    end
end
display(kontrast0);
%kontrast 45 derajat
kontrast45=0;
for i=1:ukuranx
    for j=1:ukurany
        kontrast45 = ((i-j)^2*normalisasi45(i,j))+kontrast45;
    end
end
display(kontrast45);
%kontrast 90 derajat
kontrast90=0;
for i=1:ukuranx
    for j=1:ukurany
        kontrast90 = ((i-j)^2*normalisasi90(i,j))+kontrast90;
    end
end
display(kontrast90);
%kontrast 135 derajat
kontrast135=0;
for i=1:ukuranx
    for j=1:ukurany
        kontrast135 = ((i-j)^2*normalisasi135(i,j))+kontrast135;
    end
end
display(kontrast135);
%means kontrast
kontrastm=(kontrast0+kontrast45+kontrast90+kontrast135)/4;
display(kontrastm);

%correlation 0 derajat
jml0=[0];
for i =1:ukuranx
  for j=1:ukurany
      xoi = mean(normalisasi0(i,j));
      xoj = mean(normalisasi0(i,j));
      xi = sqrt(i);
      xj = sqrt(j);
      x = (i-xoi)*(j-xoj)* (normalisasi0(i,j));
      z = x / (xi * xj);
      jml0(1,1)= jml0(1,1) + z;
  end
end

%energy 0 derajat
su0=[0];
for i =1:ukuranx
  for j=1:ukurany
  z=normalisasi0(i,j)*normalisasi0(i,j);
  su0(1,1)=su0(1,1)+z;
  end
end
energy0=su0(1,1);
display(energy0);
%energy 45 derajat
su45=[0];
for i =1:ukuranx
  for j=1:ukurany
  z=normalisasi45(i,j)*normalisasi45(i,j);
  su45(1,1)=su45(1,1)+z;
  end
end
energy45=su45(1,1);
display(energy45);
%energy 90 derajat
su90=[0];
for i =1:ukuranx
  for j=1:ukurany
  z=normalisasi90(i,j)*normalisasi90(i,j);
  su90(1,1)=su90(1,1)+z;
  end
end
energy90=su90(1,1);
display(energy90);
%energy 135 derajat
su135=[0];
for i =1:ukuranx
  for j=1:ukurany
  z=normalisasi135(i,j)*normalisasi135(i,j);
  su135(1,1)=su135(1,1)+z;
  end
end
energy135=su135(1,1);
display(energy135);
%means energy
energym=(energy0+energy45+energy90+energy135)/4;
display(energym);

correlation0=jml0(1,1);
display(correlation0);
%correlation 45 derajat
jml45=[0];
for i =1:ukuranx
  for j=1:ukurany
      xoi = mean(normalisasi45(i,j));
      xoj = mean(normalisasi45(i,j));
      xi = sqrt(i);
      xj = sqrt(j);
      x = (i-xoi)*(j-xoj)* (normalisasi45(i,j));
      z = x / (xi * xj);
      jml45(1,1)= jml45(1,1) + z;
  end
end
correlation45=jml45(1,1);
display(correlation45);
%correlation 90 derajat
jml90=[0];
for i =1:ukuranx
  for j=1:ukurany
      xoi = mean(normalisasi90(i,j));
      xoj = mean(normalisasi90(i,j));
      xi = sqrt(i);
      xj = sqrt(j);
      x = (i-xoi)*(j-xoj)* (normalisasi90(i,j));
      z = x / (xi * xj);
      jml90(1,1)= jml90(1,1) + z;
  end
end
correlation90=jml90(1,1);
display(correlation90);
%correlation 135 derajat
jml135=[0];
for i =1:ukuranx
  for j=1:ukurany
      xoi = mean(normalisasi135(i,j));
      xoj = mean(normalisasi135(i,j));
      xi = sqrt(i);
      xj = sqrt(j);
      x = (i-xoi)*(j-xoj)* (normalisasi135(i,j));
      z = x / (xi * xj);
      jml135(1,1)= jml135(1,1) + z;
  end
end
correlation135=jml135(1,1);
display(correlation135);
%means correlation
correlationm=(correlation0+correlation45+correlation90+correlation135)/4;
display(correlationm);

%homogenitas 0 derajat
hom0=[0];
for i =1:ukuranx
  for j=1:ukurany
  x=1+abs(i-j);
  z=normalisasi0(i,j)/x;
  hom0(1,1)=hom0(1,1)+z;
  end
end
homogenitas0=hom0(1,1);
display(homogenitas0);
%homogenitas 45 derajat
hom45=[0];
for i =1:ukuranx
  for j=1:ukurany
  x=1+abs(i-j);
  z=normalisasi45(i,j)/x;
  hom45(1,1)=hom45(1,1)+z;
  end
end
homogenitas45=hom45(1,1);
display(homogenitas45);
%homogenitas 90 derajat
hom90=[0];
for i =1:ukuranx
  for j=1:ukurany
  x=1+abs(i-j);
  z=normalisasi90(i,j)/x;
  hom90(1,1)=hom90(1,1)+z;
  end
end
homogenitas90=hom90(1,1);
display(homogenitas90);
%homogenitas 135 derajat
hom135=[0];
for i =1:ukuranx
  for j=1:ukurany
  x=1+abs(i-j);
  z=normalisasi135(i,j)/x;
  hom135(1,1)=hom135(1,1)+z;
  end
end
homogenitas135=hom135(1,1);
display(homogenitas135);
%means homogenitas
homogenitasm=(homogenitas0+homogenitas45+homogenitas90+homogenitas135)/4;
display(homogenitasm);

%sebagai variabel array yang akan menyimpan data training
Z1=[];
Z=[kontrast0 kontrast45 kontrast90 kontrast135 energy0 energy45 energy90 energy135 correlation0 correlation45 correlation90 correlation135 homogenitas0 homogenitas45 homogenitas90 homogenitas135];
Z1=[Z1;Z];
end
xlswrite('data_training',Z1);