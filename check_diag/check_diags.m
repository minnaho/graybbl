% requires double precision in roms his file
% change float to double in nc_read_write.opt
%gname = 'pipes_10m_hydro_impboton_grayn2_Wimp_checkdiag_grd.nc';
%fname = 'pipes_10m_hydro_impboton_grayn2_Wimp_checkdiag_his.20000101000000.nc';
%dname = 'pipes_10m_hydro_impboton_grayn2_Wimp_checkdiag_dia.20000101054500.nc';

gname = 'pipe_3m_dt1_nh_F0_debug_1ts_closed_sigma10_grd.nc';
fname = 'pipe_3m_dt1_nh_F0_debug_1ts_closed_sigma10_his.20000101000000.nc';
dname = 'pipe_3m_dt1_nh_F0_debug_1ts_closed_sigma10_dia.20000101015500.nc';

%gname = 'pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_grd.nc';
%fname = 'pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_his.20000101000000.nc';
%dname = 'pipes_10m_hydro_impboton_nogray_Wimp_checkdiag_closed_dia.20000101054500.nc';

h    = ncread(gname,'h');
msk  = ncread(gname,'mask_rho');
zeta0= ncread(fname,'zeta',[1 1 1],[inf inf 1]);
zeta1= ncread(fname,'zeta',[1 1 2],[inf inf 1]);

msku = msk(2:end,:).*msk(1:end-1,:);
msku = msku(2:end-1,2:end-1);

tim = ncread(fname,'ocean_time',[1],[2]);
delt = tim(2:end)-tim(1:end-1);
delt = sum(delt);

hc = ncreadatt(fname,'/','hc');
ts = ncreadatt(fname,'/','theta_s');
tb = ncreadatt(fname,'/','theta_b');

%N = 64;
%N = 16;
N = 10;

msku3d = repmat(msku,[1 1 N]);

zw0= permute(zlevs3(h,zeta0,ts,tb,hc,N,'w','new2008'),[2 3 1]);
zw1= permute(zlevs3(h,zeta1,ts,tb,hc,N,'w','new2008'),[2 3 1]);

dz0 = zw0(:,:,2:end)- zw0(:,:,1:end-1); 
dz1 = zw1(:,:,2:end)- zw1(:,:,1:end-1); 

% Check u diagnostics

dzu0 = 0.5*(dz0(2:end,:,:) + dz0(1:end-1,:,:));
dzu1 = 0.5*(dz1(2:end,:,:) + dz1(1:end-1,:,:));

udz0 = dzu0.*ncread(fname,'u',[1 1 1 1],[inf inf inf 1]);
udz1 = dzu1.*ncread(fname,'u',[1 1 1 2],[inf inf inf 1]);

deludz = udz1-udz0;
deludz = deludz(2:end-1,2:end-1,:).*msku3d;

 urhs1 = mean(ncread(dname,'u_pgr',[1 1 1 1],[inf inf inf 1]),4);
 urhs2 = mean(ncread(dname,'u_cor',[1 1 1 1],[inf inf inf 1]),4);
 urhs3 = mean(ncread(dname,'u_adv',[1 1 1 1],[inf inf inf 1]),4);
 urhs4 = mean(ncread(dname,'u_dis',[1 1 1 1],[inf inf inf 1]),4);
 urhs5 = mean(ncread(dname,'u_hmx',[1 1 1 1],[inf inf inf 1]),4);
 urhs6 = mean(ncread(dname,'u_vmx',[1 1 1 1],[inf inf inf 1]),4);
 urhs7 = mean(ncread(dname,'u_cpl',[1 1 1 1],[inf inf inf 1]),4);

%urhs = urhs1 + urhs2 + urhs3 + urhs4 + urhs5 + urhs6;

 urhs1 = urhs1(2:end-1,2:end-1,:).*msku3d;
 urhs2 = urhs2(2:end-1,2:end-1,:).*msku3d;
 urhs3 = urhs3(2:end-1,2:end-1,:).*msku3d;
 urhs4 = urhs4(2:end-1,2:end-1,:).*msku3d;
 urhs5 = urhs5(2:end-1,2:end-1,:).*msku3d;
 urhs6 = urhs6(2:end-1,2:end-1,:).*msku3d;
 urhs7 = urhs7(2:end-1,2:end-1,:).*msku3d;

 urhs = urhs1 + urhs2 + urhs3 + urhs4 + urhs5 + urhs6 + urhs7;

 urhs = delt*urhs;

 u_err = (deludz - urhs);
 max(abs(u_err(:)))


imagesc(squeeze(u_err(:,64,:))');axis xy;colorbar
%  imagesc(u_err(:,:,18)');axis xy;colorbar
