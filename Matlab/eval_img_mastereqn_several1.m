clear all;

cntmin=1;   % start from file no <cntmin>
cntmax=100;   % stop with file no <cntmax>
step=1;         % step through the files by <step> 
sumover=step;   % <sumover> files in the current step will be sumed up to give the next 
                % c(r,dt), so dt = step*dt(Images)
border=5;
noiselevel=0;
% size of one pixel in microns
pixel_size=200e-3;


eval_several=true;

fid=fopen('../img/4/results.txt', 'w');


Cc=50;
for nl=0:3
    resmatrix=cell(30,5);
    resmatrix{1, 1}='#images      ';
    resmatrix{1, 2}='D [um^2/s]   '; 
    resmatrix{1, 3}='dt=10ms      ';
    resmatrix{1, 4}='dt=5ms       ';
    resmatrix{1, 5}='dt=1ms       ';
    for n=0:3
        resmatrix{n*6+2,1}=20;
        resmatrix{n*6+2+1,1}=50;
        resmatrix{n*6+2+2,1}=90;
        resmatrix{n*6+2+3,1}=190;
        resmatrix{n*6+2+4,1}=500;
        resmatrix{n*6+2+5,1}=800;
        resmatrix{n*6+2,2}=50*10^(-n);
        resmatrix{n*6+2+1,2}=50*10^(-n);
        resmatrix{n*6+2+2,2}=50*10^(-n);
        resmatrix{n*6+2+3,2}=50*10^(-n);
        resmatrix{n*6+2+4,2}=50*10^(-n);
        resmatrix{n*6+2+5,2}=50*10^(-n);
    end
    cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);

    dtImages=10e-3; Dc=50;
    filename='../img/4/C50_D50_10ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn 
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{2,3}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{3,3}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{4,3}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);

    dtImages=10e-3;  Dc=5;
    filename='../img/4/C50_D5_10ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{8,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{9,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{10,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);

    dtImages=10e-3;  Dc=0.5;
    filename='../img/4/C50_D0.5_10ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    resmatrix{14,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{15,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{16,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);

    dtImages=10e-3;  Dc=0.05;
    filename='../img/4/C50_D0.05_10ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{20,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{21,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{22,4}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);

    dtImages=5e-3;  Dc=50;
    filename='../img/4/C50_D50_05ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{2,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{3,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{4,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{5,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);

    dtImages=5e-3;  Dc=5;
    filename='../img/4/C50_D5_05ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{8,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{9,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{10,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{11,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);

    dtImages=5e-3;  Dc=0.5;
    filename='../img/4/C50_D0.5_05ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{14,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{15,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{16,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=190; run eval_img_mastereqn
    resmatrix{17,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=5e-3;  Dc=0.05;
    filename='../img/4/C50_D0.05_05ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    resmatrix{20,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    resmatrix{21,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    resmatrix{22,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    resmatrix{23,5}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    
    dtImages=1e-3;  Dc=50;
    filename='../img/4/C50_D50_01ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    resmatrix{2,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    resmatrix{3,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    resmatrix{4,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    resmatrix{5,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=500; run eval_img_mastereqn
    resmatrix{6,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=800; run eval_img_mastereqn
    resmatrix{7,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);



    dtImages=1e-3;  Dc=5;
    filename='../img/4/C50_D5_01ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    resmatrix{8,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    resmatrix{9,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    resmatrix{10,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    resmatrix{11,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=500; run eval_img_mastereqn
    resmatrix{12,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=800; run eval_img_mastereqn
    resmatrix{13,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);



    dtImages=1e-3;  Dc=0.5;
    filename='../img/4/C50_D0.5_01ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    resmatrix{14,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    resmatrix{15,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    resmatrix{16,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    resmatrix{17,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=500; run eval_img_mastereqn
    resmatrix{18,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=800; run eval_img_mastereqn
    resmatrix{19,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);




    dtImages=1e-3;  Dc=0.05;
    filename='../img/4/C50_D0.05_01ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    resmatrix{20,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{21,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{22,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{23,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{24,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=500; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    resmatrix{25,6}=meanv;  cell_csv_write(sprintf('../img/4/resultmatrix_c%f_n%f.txt', Cc, noiselevel), resmatrix);
    cntmax=800; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    
    noiselevel=nl*0.1;
end

noiselevel=0;
Cc=5;

for nl=0:3
    Dc=50;
    dtImages=10e-3; 
    filename='../img/4/C5_D50_10ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn 
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=5e-3;  
    filename='../img/4/C5_D50_05ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=1e-3; 
    filename='../img/4/C5_D50_01ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=500; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=800; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);




    Dc=5;
    dtImages=10e-3; 
    filename='../img/4/C5_D5_10ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=5e-3; 
    filename='../img/4/C5_D5_05ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=1e-3; 
    filename='../img/4/C5_D5_01ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=500; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=800; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);





    Dc=0.5;
    dtImages=10e-3; 
    filename='../img/4/C5_D0.5_10ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=5e-3; 
    filename='../img/4/C5_D0.5_05ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=1e-3; 
    filename='../img/4/C5_D0.5_01ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=500; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=800; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);




    Dc=0.05;
    dtImages=10e-3; 
    filename='../img/4/C5_D0.05_10ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=5e-3; 
    filename='../img/4/C5_D0.05_05ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);

    dtImages=1e-3; 
    filename='../img/4/C5_D0.05_01ms_img_cnt%05d';
    cntmax=20; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=50; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=90; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=190; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=500; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    cntmax=800; run eval_img_mastereqn
    fprintf(fid, '"%s"; %10.3f; %10.3f; %10.3f; %10.3f; %10.5f; %10.5f\n', sprintf(filename, cntmax), Cc, dtImages, Dc,  noiselevel, meanv, stddev);
    
    noiselevel=nl*0.1;
end
fclose(fid);