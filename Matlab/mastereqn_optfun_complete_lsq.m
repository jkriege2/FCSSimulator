function y=mastereqn_optfun_complete_lsq(p, data, dt, border)

siz=size(data);
sizex=siz(1);
sizey=siz(2);

y=0;
pindex=1;
for x=border:(sizex-border)
    for y=border:(sizey-border)
        dataslice=data((x-1):(x+1),(y-1):(y+1),:);
        y=y+mastereqn_optfun_lsq(p(pindex), dataslice, dt);
        pindex=pindex+1;
    end
end