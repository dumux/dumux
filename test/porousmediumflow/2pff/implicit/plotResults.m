clear all;

pathToPVD = '';
testName = 'test2_constvelff_nohybrid';
pathAndName = strcat(pathToPVD,testName);

%fileIndex = {};
%lastIndex = 10;
%for k=0:lastIndex
%    fileIndex{end+1} = sprintf('%05s', num2str(k))
%end

fileIndex = {'00008','00031','00166'};

numFliles = max(size(fileIndex));
for i=1:numFliles
    fileName = strcat(pathAndName,'-',fileIndex{i},'.vtu');
    system(['pvpython ../../../../../bin/postprocessing/extractlinedata.py -f', ' ',fileName , ' ', '-p1 0.5 100 0  -p2 0.5 0 0'])
    outName = strcat(pathAndName,'-',fileIndex{i},'.csv');
    data = csvread(outName,2.0);
    Sw = data(:,1);
    Sn = data(:,2);
    z = data(:,14);

    figure;
    plot(z,Sn,'b','LineWidth',3 );
    xlabel('$z$','Interpreter','latex');
    ylabel('$S_n$','Interpreter','latex');
    ax = gca;
    grid  on;
    set(gca,'FontSize', 18);
    set(findall(gcf,'type','text'),'fontSize',28); 
    
    fnam=[testName '_Sn_' fileIndex{i} '.pdf'];
    saveas(gcf,fnam,'pdf')  
    system (['/usr/bin/pdfcrop -margin 10', ' ', fnam, ' ', fnam]); 
end
