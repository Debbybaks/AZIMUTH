// AU, day, solar mass internally. All interface conversions explicitly labelled.
export const AU=149597870700, DAY=86400, MU=0.0002959122082855911, TAU=Math.PI*2;
export const add=(a,b)=>a.map((v,i)=>v+b[i]), sub=(a,b)=>a.map((v,i)=>v-b[i]), mul=(a,s)=>a.map(v=>v*s), dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0), norm=a=>Math.hypot(...a), cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
export const rad=x=>x*Math.PI/180, deg=x=>x*180/Math.PI;
export const bodies=[
 {name:'Sun',a:0,radius:696340,mu:1.32712440018e20,color:'#ffbd61',mass:1,period:0},
 {name:'Mercury',a:.387099,e:.205636,i:7.005,node:48.339,peri:77.458,L:252.251,radius:2439.7,color:'#b7a99b',mass:1.6601e-7},
 {name:'Venus',a:.723336,e:.006777,i:3.395,node:76.68,peri:131.602,L:181.98,radius:6051.8,color:'#dcbe8a',mass:2.4478e-6},
 {name:'Earth',a:1.000003,e:.016711,i:0,node:0,peri:102.938,L:100.464,radius:6371,color:'#6cafff',mass:3.0035e-6},
 {name:'Moon',a:.00257,e:.0549,i:5.145,node:125.08,peri:83.353,L:218.316,radius:1737.4,color:'#c0c5cd',mass:3.694e-8,parent:'Earth'},
 {name:'Mars',a:1.52371,e:.093394,i:1.85,node:49.559,peri:336.041,L:355.453,radius:3389.5,color:'#dc8769',mass:3.2272e-7},
 {name:'Jupiter',a:5.202887,e:.048386,i:1.304,node:100.474,peri:14.728,L:34.396,radius:69911,color:'#c9af91',mass:.0009543},
 {name:'Saturn',a:9.536676,e:.053862,i:2.486,node:113.662,peri:92.599,L:49.954,radius:58232,color:'#d1bd8e',mass:.0002857},
 {name:'Uranus',a:19.189165,e:.047257,i:.773,node:74.016,peri:170.954,L:313.238,radius:25362,color:'#9acfd3',mass:.00004365},
 {name:'Neptune',a:30.069923,e:.00859,i:1.77,node:131.784,peri:44.965,L:304.88,radius:24622,color:'#5789e2',mass:.00005149}
].map(b=>({...b,period:b.a?TAU*Math.sqrt(b.a**3/(MU*(b.parent?3.0035e-6:1))):0}));
export const epochDays=date=>(new Date(date).getTime()-Date.UTC(2000,0,1,12))/86400000;
export function elementsState({a,e,i=0,node=0,omega=0,nu=0},mu=MU){
 if(!(a>0&&e>=0&&e<1))throw Error('Elliptic orbit requires a > 0 and 0 ≤ e < 1.');
 const p=a*(1-e*e),r=p/(1+e*Math.cos(nu));
 const c=Math.cos,s=Math.sin,O=node,w=omega,I=i;
 const P=[c(O)*c(w)-s(O)*s(w)*c(I),s(O)*c(w)+c(O)*s(w)*c(I),s(w)*s(I)];
 const Q=[-c(O)*s(w)-s(O)*c(w)*c(I),-s(O)*s(w)+c(O)*c(w)*c(I),c(w)*s(I)];
 return {r:add(mul(P,r*c(nu)),mul(Q,r*s(nu))),v:add(mul(P,-Math.sqrt(mu/p)*s(nu)),mul(Q,Math.sqrt(mu/p)*(e+c(nu))))};
}
export function bodyState(name,t){
 const b=bodies.find(b=>b.name===name);if(!b)throw Error('Unknown body');if(!b.a)return {r:[0,0,0],v:[0,0,0]};
 const mu=MU*(b.parent?3.0035e-6:1),M=rad(b.L-b.peri)+Math.sqrt(mu/b.a**3)*t;
 let E=M;for(let j=0;j<15;j++)E-=(E-b.e*Math.sin(E)-M)/(1-b.e*Math.cos(E));
 const nu=2*Math.atan2(Math.sqrt(1+b.e)*Math.sin(E/2),Math.sqrt(1-b.e)*Math.cos(E/2));
 const state=elementsState({a:b.a,e:b.e,i:rad(b.i),node:rad(b.node),omega:rad(b.peri-b.node),nu},mu);
 if(b.parent){const parent=bodyState(b.parent,t);state.r=add(state.r,parent.r);state.v=add(state.v,parent.v);}return state;
}
export function stateElements(r,v,mu=MU){
 const R=norm(r),V=norm(v),h=cross(r,v),H=norm(h),n=cross([0,0,1],h),N=norm(n),ev=sub(mul(cross(v,h),1/mu),mul(r,1/R)),e=norm(ev),energy=V*V/2-mu/R,a=-mu/(2*energy),clamp=x=>Math.max(-1,Math.min(1,x));
 const i=Math.acos(clamp(h[2]/H)),node=N?Math.atan2(n[1],n[0]):0;
 let omega=N&&e>1e-10?Math.acos(clamp(dot(n,ev)/(N*e))):Math.atan2(ev[1],ev[0]);if(N&&ev[2]<0)omega=TAU-omega;
 let nu=e>1e-10?Math.acos(clamp(dot(ev,r)/(e*R))):Math.atan2(r[1],r[0]);if(e>1e-10&&dot(r,v)<0)nu=TAU-nu;
 return {a,e,i,node,omega,nu,energy,h:H,period:a>0?TAU*Math.sqrt(a**3/mu):null,periapsis:a*(1-e),apoapsis:e<1?a*(1+e):null};
}
export function hohmann(r1,r2,mu=MU){if(!(r1>0&&r2>0))throw Error('Positive radii required');const a=(r1+r2)/2,dv1=Math.sqrt(mu/r1)*(Math.sqrt(2*r2/(r1+r2))-1),dv2=Math.sqrt(mu/r2)*(1-Math.sqrt(2*r1/(r1+r2)));return {dv1,dv2,total:Math.abs(dv1)+Math.abs(dv2),time:Math.PI*Math.sqrt(a**3/mu),a,e:Math.abs(r2-r1)/(r1+r2)};}
function stumpff(z){if(z>1e-7){const q=Math.sqrt(z);return [(1-Math.cos(q))/z,(q-Math.sin(q))/q**3];}if(z< -1e-7){const q=Math.sqrt(-z);return [(Math.cosh(q)-1)/(-z),(Math.sinh(q)-q)/q**3];}return [.5-z/24+z*z/720,1/6-z/120+z*z/5040];}
export function lambert(r1,r2,tof,mu=MU){
 if(!(tof>0))throw Error('Time of flight must be positive');const R1=norm(r1),R2=norm(r2),cos=dot(r1,r2)/(R1*R2),cz=cross(r1,r2)[2],sin=Math.sign(cz||1)*Math.sqrt(Math.max(0,1-cos*cos));const A=sin*Math.sqrt(R1*R2/(1-cos));
 if(!Number.isFinite(A)||Math.abs(A)<1e-8)throw Error('Lambert geometry is singular near 0° or 180°. Adjust epoch or duration.');
 const f=z=>{const [C,S]=stumpff(z);if(C<=0)return null;const y=R1+R2+A*(z*S-1)/Math.sqrt(C);if(y<0)return null;return {y,t:((y/C)**1.5*S+A*Math.sqrt(y))/Math.sqrt(mu)};};
 let lo=-4*Math.PI*Math.PI,hi=4*Math.PI*Math.PI-1e-5,z=0,q,iterations=0;
 for(;iterations<160;iterations++){z=(lo+hi)/2;q=f(z);if(!q||q.t<tof)lo=z;else hi=z;if(q&&Math.abs(q.t-tof)<1e-8)break;}
 if(!q||Math.abs(q.t-tof)>1e-5)throw Error('Lambert did not converge within 160 iterations.');
 const g=A*Math.sqrt(q.y/mu),F=1-q.y/R1,G=1-q.y/R2;
 return {v1:mul(sub(r2,mul(r1,F)),1/g),v2:mul(sub(mul(r2,G),r1),1/g),iterations,residual:Math.abs(q.t-tof),angle:Math.atan2(sin,cos)};
}
export const acceleration=(r,mu=MU)=>mul(r,-mu/norm(r)**3);
export function step(r,v,h,method='RK4',mu=MU,force=null){
 const acc=force||((r)=>acceleration(r,mu)),a=acc(r);
 if(method==='Euler')return {r:add(r,mul(v,h)),v:add(v,mul(a,h)),evals:1};
 if(method==='Verlet'){const rn=add(add(r,mul(v,h)),mul(a,h*h/2));return {r:rn,v:add(v,mul(add(a,acc(rn)),h/2)),evals:2};}
 if(method==='Heun'){const rp=add(r,mul(v,h)),vp=add(v,mul(a,h));return {r:add(r,mul(add(v,vp),h/2)),v:add(v,mul(add(a,acc(rp)),h/2)),evals:2};}
 const k2r=add(v,mul(a,h/2)),k2v=acc(add(r,mul(v,h/2))),k3r=add(v,mul(k2v,h/2)),k3v=acc(add(r,mul(k2r,h/2))),k4r=add(v,mul(k3v,h)),k4v=acc(add(r,mul(k3r,h)));
 return {r:add(r,mul(add(add(v,mul(k2r,2)),add(mul(k3r,2),k4r)),h/6)),v:add(v,mul(add(add(a,mul(k2v,2)),add(mul(k3v,2),k4v)),h/6)),evals:4};
}
export function missionAcceleration(pos,epoch,model='Two-body'){let a=acceleration(pos);if(model==='Restricted N-body')for(const b of bodies.filter(b=>b.a&&!b.parent)){const br=bodyState(b.name,epoch).r,d=sub(br,pos);if(norm(d)>b.radius*1000/AU*2)a=add(a,sub(mul(d,MU*b.mass/norm(d)**3),mul(br,MU*b.mass/norm(br)**3)));}return a;}
export function simulate(config,override=null){
 const {departure,destination,date,duration,solver='RK4',stepSize=.5,burn=0,angle=0,model='Two-body'}=config;
 if(departure===destination)throw Error('Choose different departure and destination bodies.');
 if(!(duration>0&&duration<=2000&&stepSize>=.01&&stepSize<=10))throw Error('Duration 1–2000 days; integration step 0.01–10 days.');
 const t0=epochDays(date),start=bodyState(departure,t0),end=bodyState(destination,t0+duration),solution=lambert(start.r,end.r,duration);
 let initial=override||add(solution.v1,[burn*DAY/AU*Math.cos(rad(angle)),burn*DAY/AU*Math.sin(rad(angle)),0]);
 let r=[...start.r],v=[...initial],t=0,points=[],closest=Infinity,closestTime=0,evals=0,maxDrift=0;const energy0=norm(v)**2/2-MU/norm(r),beg=performance.now();
 const force=model==='Restricted N-body'?(pos)=>missionAcceleration(pos,t0+t,model):null;
 const maxSteps=Math.ceil(duration/stepSize);if(maxSteps>200000)throw Error('Step budget exceeded');
 for(let j=0;j<=maxSteps;j++){const target=bodyState(destination,t0+t),distance=norm(sub(r,target.r));if(distance<closest){closest=distance;closestTime=t;}const E=norm(v)**2/2-MU/norm(r),drift=Math.abs((E-energy0)/energy0);maxDrift=Math.max(maxDrift,drift);
 if(j%Math.max(1,Math.floor(maxSteps/1600))===0||t>=duration)points.push({t,r:[...r],v:[...v],distance,energy:E,drift});
 if(t>=duration)break;const h=Math.min(stepSize,duration-t);const next=step(r,v,h,solver,MU,force);r=next.r;v=next.v;t+=h;evals+=next.evals;if(!Number.isFinite(norm(r))||norm(r)<.001)throw Error('Numerical failure or solar collision. Reduce the step size.');}
 const miss=norm(sub(r,end.r)),arrival=norm(sub(v,end.v)),dv=norm(sub(initial,start.v));
 return {config,points,initial,solution,closest,closestTime,miss,arrival,dv,energy0,maxDrift,evals,computeMs:performance.now()-beg,status:miss*AU<1e8?'Encounter':'Miss',events:[{t:0,name:'Departure burn'},{t:closestTime,name:'Closest approach'},{t:duration,name:miss*AU<1e8?'Encounter (100,000 km tolerance)':'Arrival / miss'}],assumptions:'Fixed J2000 Kepler elements; heliocentric inertial ecliptic frame. Impulsive departure, no parking-orbit escape or arrival capture. Lambert uses solar two-body gravity; N-body propagation adds approximate planetary perturbations. No operational ephemerides.'};
}
export function coordinates(x,y,z,lat=0,lst=0){const r=Math.hypot(x,y,z);if(!r)throw Error('The zero vector has no angular direction');const ra=Math.atan2(y,x),dec=Math.asin(z/r),eps=rad(23.43928),eq=[x,y*Math.cos(eps)-z*Math.sin(eps),y*Math.sin(eps)+z*Math.cos(eps)],R=Math.atan2(eq[1],eq[0]),D=Math.asin(eq[2]/r),H=rad(lst)-R,p=rad(lat),alt=Math.asin(Math.sin(p)*Math.sin(D)+Math.cos(p)*Math.cos(D)*Math.cos(H)),az=Math.atan2(-Math.sin(H)*Math.cos(D),Math.sin(D)*Math.cos(p)-Math.cos(D)*Math.sin(p)*Math.cos(H));return {r,longitude:(deg(ra)+360)%360,latitude:deg(dec),ra:((deg(R)+360)%360)/15,dec:deg(D),alt:deg(alt),az:(deg(az)+360)%360,equatorial:eq};}
export function gravityAssist(vinf,periapsis,mu=3.986004418e14,planetSpeed=29780,heading=0){if(vinf<=0||periapsis<6371000)throw Error('Positive excess speed and periapsis above Earth surface required');const e=1+periapsis*vinf**2/mu,bend=2*Math.asin(1/e),h=rad(heading),incoming=[vinf*Math.cos(h)+planetSpeed,vinf*Math.sin(h),0],outgoing=[vinf*Math.cos(h+bend)+planetSpeed,vinf*Math.sin(h+bend),0];return {e,bend,initial:norm(incoming),final:norm(outgoing),change:norm(outgoing)-norm(incoming)};}
export function groundTrack(altitude=550,inclination=53,count=1,spacing=0){const radius=6371000+altitude*1000,mu=3.986004418e14,period=TAU*Math.sqrt(radius**3/mu),incl=rad(inclination);let tracks=[];for(let k=0;k<count;k++){let points=[];for(let j=0;j<=360;j++){const t=j/360*period*2,nu=TAU*t/period+k*TAU/count,lon=Math.atan2(Math.sin(nu)*Math.cos(incl),Math.cos(nu))-TAU*t/86164+rad(spacing*k),lat=Math.asin(Math.sin(nu)*Math.sin(incl));points.push({lon:((deg(lon)+540)%360)-180,lat:deg(lat),t});}tracks.push(points);}return {period,tracks,footprint:Math.acos(6371000/radius),radius};}
export function decay(altitude,mass,area,cd=2.2){const h=altitude*1000,r=6371000+h,rho=3.614e-13*Math.exp(-(h-700000)/88667),v=Math.sqrt(3.986004418e14/r),a=.5*rho*cd*area/mass*v*v,rate=-2*a*r/v;return {rho,drag:a,rate:rate*DAY,period:TAU*r/v};}
// Fully interacting planar three-body experiment. Unlike restricted mission mode,
// these states evolve together under mutual Newtonian acceleration.
export function nbody({duration=800,stepSize=1,massMultiplier=1}={}){
 if(!(duration>=10&&duration<=3000&&stepSize>=.05&&stepSize<=5&&massMultiplier>=.1&&massMultiplier<=20))throw Error('N-body limits: 10–3000 days, 0.05–5 day steps, mass multiplier 0.1–20.');
 const masses=[1,.0009543*massMultiplier,3.0035e-6],r=[[0,0,0],[5.2,0,0],[1,0,0]],v=[[0,0,0],[0,Math.sqrt(MU/5.2),0],[0,Math.sqrt(MU),0]],total=masses.reduce((a,b)=>a+b,0);
 const com=r.reduce((a,x,i)=>add(a,mul(x,masses[i]/total)),[0,0,0]),momentum=v.reduce((a,x,i)=>add(a,mul(x,masses[i]/total)),[0,0,0]);for(let i=0;i<3;i++){r[i]=sub(r[i],com);v[i]=sub(v[i],momentum);}
 const force=()=>r.map((ri,i)=>r.reduce((a,rj,j)=>{if(i===j)return a;const d=sub(rj,ri),R=norm(d);if(R<.0001)throw Error('Close collision: integration stopped.');return add(a,mul(d,MU*masses[j]/R**3));},[0,0,0]));
 const energy=()=>{let e=0;for(let i=0;i<3;i++){e+=masses[i]*norm(v[i])**2/2;for(let j=i+1;j<3;j++)e-=MU*masses[i]*masses[j]/norm(sub(r[i],r[j]));}return e;};
 const e0=energy(),points=[];let a=force(),t=0,maxDrift=0;const stride=Math.max(1,Math.ceil(duration/stepSize/1200));let j=0;
 while(t<duration){const h=Math.min(stepSize,duration-t);for(let i=0;i<3;i++)r[i]=add(add(r[i],mul(v[i],h)),mul(a[i],h*h/2));const an=force();for(let i=0;i<3;i++)v[i]=add(v[i],mul(add(a[i],an[i]),h/2));a=an;t+=h;const drift=Math.abs((energy()-e0)/e0);maxDrift=Math.max(drift,maxDrift);if(j++%stride===0||t>=duration)points.push({t,r:r.map(x=>[...x]),drift});if(r.some(x=>!Number.isFinite(norm(x))||norm(x)>1000))throw Error('Runaway state: experiment stopped.');}
 return {points,masses,maxDrift,steps:j,method:'Velocity Verlet',units:'AU, day, solar mass'};
}
export function solverArena({stepSize=4,orbits=3}={}){
 if(!(stepSize>=.05&&stepSize<=50&&orbits>=1&&orbits<=20))throw Error('Use 0.05–50 day steps and 1–20 orbits');
 const period=TAU/Math.sqrt(MU),duration=period*orbits,stride=Math.max(1,Math.ceil(duration/stepSize/1200)),rows=[],series=[];
 for(const [method,color]of [['Euler','#e5a471'],['Heun','#75a9e8'],['RK4','#62dfc0'],['Verlet','#ad91ed']]){let r=[1,0,0],v=[0,Math.sqrt(MU),0],t=0,n=0,j=0,points=[],max=0;const start=performance.now();while(t<duration){const h=Math.min(stepSize,duration-t),s=step(r,v,h,method);r=s.r;v=s.v;t+=h;n+=s.evals;const drift=Math.abs((norm(v)**2/2-MU/norm(r)+MU/2)/(MU/2));max=Math.max(max,drift);if(j++%stride===0||t>=duration)points.push([t,Math.log10(Math.max(drift,1e-16))]);}const expected=[Math.cos(Math.sqrt(MU)*duration),Math.sin(Math.sqrt(MU)*duration),0];rows.push({method,error:norm(sub(r,expected)),drift:max,evals:n,ms:performance.now()-start});series.push({points,color,label:method});}return {rows,series};
}
