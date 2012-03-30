function littlerock_ode45
   %plot a sounding 
    filename='littlerock.nc';
    fprintf('reading file: %s\n',filename);
    file_struct=nc_info(filename);
    c=constants;
    %
    % grab the March 2 12Z sounding
    %
    sound_var = file_struct.Dataset(4).Name;
    fprintf('found sounding: %s\n',sound_var);
    press=nc_varget(filename,sound_var,[0,0],[Inf,1]);
    temp=nc_varget(filename,sound_var,[0,2],[Inf,1]);
    dewpoint=nc_varget(filename,sound_var,[0,3],[Inf,1]);
    height=nc_varget(filename,sound_var,[0,1],[Inf,1]);
    
    newHeight=nudgeHeight(height);
     
    InterpTenv=@(zVals) interp1(newHeight,temp,zVals);
    InterpTdEnv=@(zVals) interp1(newHeight, dewpoint,zVals);
    InterpPress=@(zVals) interp1(newHeight,press,zVals);
   
    P900=find(abs(900-press) < 2.);
    P800=find(abs(800-press) < 7.);
    height_800=height(P800);
    thetae=thetaep(dewpoint(P900) + c.Tc,temp(P900)+c.Tc,press(P900)*100.);
       
    tspan=0:10:2500;
    yinit=[0.5,height_800];
    derivs=@(t,y) F(t,y,thetae,InterpTenv,InterpTdEnv,InterpPress);
    [t,y]=ode45(derivs,tspan,yinit);
    wvel=y(:,1);
    height=y(:,2);
    
    updraft=wvel>0;
    figure(1);
    clf;
    plot(wvel(updraft),height(updraft))
    xlabel('vertical velocity');
    ylabel('height above surface');
    title('wvel vs. height')
    grid on;
end   

  function yp=F(t,y,thetae0,InterpTenv,InterpTdEnv,InterpPress)
  yp=zeros(2,1); % since output must be a column vector
  yp(1)=calcBuoy(y(2),thetae0,InterpTenv,InterpTdEnv,InterpPress); 
  yp(2)= y(1);
  end 
  
function Bout=calcBuoy(height,thetae0,InterpTenv,InterpTdEnv,InterpPress)
    %calcTvDiff(press,thetae0,interpTenv,interpTdenv)
    %input: press (Pa), thetae0 (K), plus function handles for T,Td soundings
    %output: TvDiff (K)
    %neglect liquid water loading in the virtual temperature
    c=constants;
    press=InterpPress(height)*100;
    Tcloud=findTmoist(thetae0,press);
    wvcloud=wsat(Tcloud,press);
    Tvcloud=Tcloud*(1. + c.eps*wvcloud);
    Tenv=InterpTenv(press*1.e-2) + c.Tc;
    Tdenv=InterpTdEnv(press*1.e-2) + c.Tc;
    wvenv=wsat(Tdenv,press);
    Tvenv=Tenv*(1. + c.eps*wvenv);
    TvDiff=Tvcloud - Tvenv;
    fprintf('%10.3f %10.3f %10.3f\n',press*0.01,height,TvDiff);
    Bout=c.g0*(TvDiff/Tvenv);
end



function newHeight=nudgeHeight(zVec)
     %if two balloon height levels are idential or descend
     %add a factor of 0.1% to the second one
     %so interpolation will work
     newHeight=zVec;
     hit=find(diff(newHeight) <= 0);
     newHeight(hit+1)=zVec(hit) + 1.e-3*zVec(hit);
end

