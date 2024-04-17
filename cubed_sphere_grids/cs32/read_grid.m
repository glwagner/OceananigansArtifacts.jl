
nc=32; np=nc+1; nx=6*nc; nPg=nx*nc;
if size(who('krd'),1) > 0,
 fprintf('krd is defined and = %i \n',krd);
else
 fprintf('krd undefined ; set to 1 \n'); krd=1 ;
end

inpName=['grid_cs',int2str(nc)];
%rDir='/home/jmc/exp/inp_CS_grid/';
rDir='./';

if krd > 0,
 for n=1:6,
%-- read :
  namF=sprintf('%s.%s%3.3i.%s',[rDir,inpName],'face',n,'bin');
  fid=fopen(namF,'r','b');
  vv1=fread(fid,'real*8');
  fclose(fid);
  s=size(vv1,1);
  fprintf(['read: ',namF,' : size: %i (%ix%ix%i)\n'],s,np,np,s/np/np);
  k=s/np/np;
  vv1=reshape(vv1,[np np k]);
  if n == 1, vv6=zeros(np,np,k,6); end
  vv6(:,:,:,n)=vv1;
 end
end

titv={'xC','yC','dxF','dyF','rAc','xG','yG','dxV','dyU','rAz','dxC','dyC', ...
      'rAw','rAs','dxG','dyG','AngleCS','AngleSN'};
%- XC YC DXF DYF RA XG YG DXV DYU RAZ DXC DYC RAW RAS DXG DYG AngleCS AngleSN :
%------------------------------------------------------------------------------

%- make a plot:
xax=[1:np]-0.5; yax=xax;
%nv=1;  ccB=[-1 1]*180;
%nv=2;  ccB=[-1 1]*90;
%nv=3;  ccB=[1 3.2]*1.e+5;
%nv=4;  ccB=[1 3.2]*1.e+5;
%nv=5;  ccB=[1 10.]*1.e+10;
%nv=6;  ccB=[-1 1]*180;
 nv=7;  ccB=[-1 1]*90;
%nv=8;  ccB=[1 3.2]*1.e+5;
%nv=9;  ccB=[1 3.2]*1.e+5;
%nv=10; ccB=[1 10.]*1.e+10;
%nv=11; ccB=[1 3.2]*1.e+5;
%nv=12; ccB=[1 3.2]*1.e+5;
%nv=13; ccB=[1 10.]*1.e+10;
%nv=14; ccB=[1 10.]*1.e+10;
%nv=15; ccB=[1 3.2]*1.e+5;
%nv=16; ccB=[1 3.2]*1.e+5;
%nv=17; ccB=[-1 1];
%nv=18; ccB=[-1 1];

%- default: dxRed=0; dyRed=0; dxB=0.1; dyB=0.9;
dxRed=0.05; dyRed=0; dxB=0.02; dyB=0.9;
[xyP,xyB]=def_subP(6,dxRed,dyRed,dxB,dyB);

figure(1); clf; colormap('jet');
for j=1:6,
 var=squeeze(vv6(:,:,nv,j));
 axes('position',xyP(j,:));
 imagesc(xax,yax,var'); set(gca,'YDir','normal');
 caxis(ccB);
 if j==3,
   BB=colorbar; set(BB,'Position',xyB(j,:));
 end
 title([char(titv(nv)),' , face: ',int2str(j)]);
end

