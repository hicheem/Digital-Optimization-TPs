set of int: v_rnd={3,6,1,8,4,9,0,5,7,2};
%output ["vrand = \(v_rnd);"];

array[0..9] of var v_rnd: v;
array[0..4] of var v_rnd: vp;
array[0..4] of var v_rnd: vimp;

%output ["v = \(v);"];
constraint forall(e in v_rnd)(exists(i in 0..9)(v[i]=e));
constraint forall(i,j in 0..9 where i<j )(v[i] <= v[j]);
constraint forall(e in v where e mod 2==0)(exists(k in 0..4)(vp[k] = e));
constraint forall(e in v where e mod 2!=0)(exists(k in 0..4)(vimp[k] = e));
constraint forall(i,j in 0..4 where i<j )(vp[i] <= vp[j]);
constraint forall(i,j in 0..4 where i<j )(vimp[i] <= vimp[j]);
solve satisfy;

output[ "v=\(v);\n"];
output["v_impaire = \(vimp)\n"++"v_paire=\(vp)"];