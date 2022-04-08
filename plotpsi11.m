testDir='F:\cbwu\个人信息\中文核心\newdata for cfetr 20181206';
%gfile = fullfile(testDir,'g080307.040000_old');
gfile = fullfile(testDir,'g888888.888888');
gfileVar = readGfile1(gfile); 
gvar = gfileVar;

%%
X = linspace(0, gvar.rdim, gvar.nw);
Z = linspace(0, gvar.zdim, gvar.nw);
rgrid = X+gvar.rleft;
zgrid = Z-(gvar.zmid+gvar.zdim/2);
psirz = gvar.psirz;
ssimag = gvar.simag;    %磁轴处psi值
ssibry = gvar.sibry;    %    最外闭合磁面处psi值
pres = gvar.pres;
%%
a=(max(gvar.rbbbs)-min(gvar.rbbbs))/2;
b=(max(gvar.zbbbs)-min(gvar.zbbbs))/2;
R=(max(gvar.rbbbs)+min(gvar.rbbbs))/2;
k=b/a;
dtriu=(R-gvar.rbbbs(find(gvar.zbbbs==max(gvar.zbbbs))))/a;
dtrid=(R-gvar.rbbbs(find(gvar.zbbbs==min(gvar.zbbbs))))/a;
rhopsi=(psirz'-ssimag)/(ssibry-ssimag);

% rhopsip=interp2(rgrid,zgrid,rhopsi,rgrid,zgrid,'cubic');

 
%%
figure;
%contour(rgrid,zgrid,-gvar.psirz',20,'r','LineWidth',1.5);
contourf(rgrid,zgrid,-gvar.psirz',50)

hold on;
plot(gvar.rbbbs,gvar.zbbbs,'r','LineWidth',2);   %最外闭合磁面
plot(gvar.rlim,gvar.zlim,'b','LineWidth',2);    %限制器
% colorbar();
axis equal;
xlabel('R(m)');
ylabel('Z(m)');
hold on;

function gVar = readGfile1( fileName, replaceMat)
%Read gfile by jlchen.
% 
% fileName = 'g034069.09999'
% whos fileName
 %#ok<*NASGU>
if nargin<2
    replaceMat = 0;
end 

if ~exist(fileName, 'file')
    gVar = [];
    warning(['gfile:',  fileName, ' is not found']);
    return;
end
matFileName = [fileName,'_tmp_', '.mat'];
if ~exist(matFileName, 'file') || replaceMat ~= 0     
    fid = fopen(fileName);
    firstLineDFormat = '%4d'; 
    eFormat = '%e'; eNum = 5;
    dFormat = '%5d';
    case6 = transpose(fscanf(fid, '%c', [8,6]));
    % C = textscan(fid, '%s Level%u8 %f32 %d8 %u %f %f %s %f', 1);
    [idum, nw, nh] = fscanf2each(fid, '%d',3); %#ok<*ASGLU>
    % idum
    % nw
    % nh
    % fscanf(fid, '%e', eNum)
    % fscanf(fid, eFormat, eNum)
    fgetl(fid);
    [rdim,zdim,rcentr,rleft,zmid] = fscanf2each(fid, eFormat, eNum);
    % rdim,zdim,rcentr,rleft,zmid
    [rmaxis,zmaxis,simag,sibry,bcentr] = fscanf2each(fid, eFormat, eNum);
    [current,simag,xdum,rmaxis,xdum] = fscanf2each(fid, eFormat, eNum);
    [zmaxis,xdum,sibry,xdum,xdum] = fscanf2each(fid, eFormat, eNum);
    fpol = fscanf(fid, eFormat, nw);
    pres = fscanf(fid, eFormat, nw);
    ffprime = fscanf(fid, eFormat, nw);
    pprime = fscanf(fid, eFormat, nw);
    psirz = fscanf(fid, eFormat, [nw, nh]);
    qpsi = fscanf(fid, eFormat, nw);

    nbbbs = fscanf(fid, '%d', 1);
    limitr  = fscanf(fid, '%d',1);
    tmp = fscanf(fid, eFormat, [2, nbbbs]);
    rbbbs = tmp(1,:);
    zbbbs = tmp(2,:);
    tmp = fscanf(fid, eFormat, [2, limitr]);
    rlim = tmp(1,:);
    zlim = tmp(2,:);

    % character*10 case(6)
    % dimension psirz(nw,nh),fpol(1),pres(1),ffprim(1),
    % . pprime(1),qpsi(1),rbbbs(1),zbbbs(1),
    % . rlim(1),zlim(1)
    % c
    % read (neqdsk,2000) (case(i),i=1,6),idum,nw,nh
    % read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
    % read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
    % read (neqdsk,2020) current,simag,xdum,rmaxis,xdum
    % read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
    % read (neqdsk,2020) (fpol(i),i=1,nw)
    % read (neqdsk,2020) (pres(i),i=1,nw)
    % read (neqdsk,2020) (ffprim(i),i=1,nw)
    % read (neqdsk,2020) (pprime(i),i=1,nw)
    % read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
    % read (neqdsk,2020) (qpsi(i),i=1,nw)
    % read (neqdsk,2022) nbbbs,limitr
    % read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
    % read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
    % c
    % 2000 format (6a8,3i4)
    % 2020 format (5e16.9)
    % 2022 format (2i5)
    fclose(fid);

    clear eFormat dFormat firstLineDFormat 
    clear tmp
    clear fid
    
    save(matFileName);
%     gVar = load(['gfiel_tmp_', fileName, '.mat']);
%     delete(['gfiel_tmp_', fileName, '.mat']);
end
gVar = load(matFileName);

end

function varargout = fscanf2each(fid, format, num)
    if nargout == 1;
        varargout = fscanf(fid, format, num);
    else
        varargout = cell(1, nargout);
        output1 = fscanf(fid, format, num);
        for jj = 1:length(output1)
            varargout{jj} = output1(jj);
        end
    end
end




