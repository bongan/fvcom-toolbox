[x,y,nx,ny,val] = textread('vec_check.dat','%f %f %f %f %f')
[bx,by,ii] = textread('jordan_basin_bndr.dat','%f %f %d')

quiver(x,y,nx,ny)
hold on
plot(bx,by,'g-+')
scatter(x,y,50,val)
