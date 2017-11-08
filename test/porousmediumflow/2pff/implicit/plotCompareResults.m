clear all;

pathToPVD = '';
testName1 = 'test2_constvelff_hybrid';
pathAndName1 = strcat(pathToPVD,testName1);

testName2 = 'test2_constvelff_nohybrid';
pathAndName2 = strcat(pathToPVD,testName2);

%fileIndex = {};
%lastIndex = 10;
%for k=0:lastIndex
%    fileIndex{end+1} = sprintf('%05s', num2str(k))
%end

fileIndex = {'00008','00031','00166'};

numFliles = max(size(fileIndex));
for i=1:numFliles
    fileName = strcat(pathAndName1,'-',fileIndex{i},'.vtu');
    system(['pvpython ../../../../../bin/postprocessing/extractlinedata.py -f', ' ',fileName , ' ', '-p1 0.5 100 0  -p2 0.5 0 0'])
    outName = strcat(pathAndName1,'-',fileIndex{i},'.csv');
    data = csvread(outName,2.0);
    Sw1 = data(:,1);
    Sn1 = data(:,2);
    z1 = data(:,14);

    fileName = strcat(pathAndName2,'-',fileIndex{i},'.vtu');
    system(['pvpython ../../../../../bin/postprocessing/extractlinedata.py -f', ' ',fileName , ' ', '-p1 0.5 100 0  -p2 0.5 0 0'])
    outName = strcat(pathAndName2,'-',fileIndex{i},'.csv');
    data = csvread(outName,2.0);
    Sw2 = data(:,1);
    Sn2 = data(:,2);
    z2 = data(:,14);

    figure;
    plot(z1,Sn1,'b','LineWidth',3 );
    hold on;
    plot(z2,Sn2,'r','LineWidth',3 );
    xlabel('$z$','Interpreter','latex');
    ylabel('$S_n$','Interpreter','latex');
    ax = gca;
    grid  on;
    set(gca,'FontSize', 18);
    set(findall(gcf,'type','text'),'fontSize',28); 
    
    fnam=[testName1 'Compare_Sn_' fileIndex{i} '.pdf'];
    saveas(gcf,fnam,'pdf')  
    system (['/usr/bin/pdfcrop -margin 10', ' ', fnam, ' ', fnam]); 
end
