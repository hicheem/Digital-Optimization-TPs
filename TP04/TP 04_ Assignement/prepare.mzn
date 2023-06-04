int: n;
set of int: TASK = 1..n;
array[TASK] of int:  d;  % length of task

int: m;
set of int: PREC = 1..m;
array[PREC] of TASK:  pre;  % task pre[i] must finish before post[t] starts
array[PREC] of TASK:  post; % 

int: o;  % This is not a 0, this the letter "o" 
set of int: DISJ = 1..o;
array[DISJ] of set of TASK: disj; % sets of tasks which can not overlap in time

int: maxt = sum(t in TASK)(d[t]);

array[TASK] of var 0..maxt: s; % start time
array[TASK] of var 0..maxt: e = [ s[t] + d[t] | t in TASK]; % end time


% performs tasks obey precedences
constraint
   forall ( i in PREC, j in pre )
   (
         s[i] + d[i] >= s[j]
   );

% no two disjoint tasks overlap
constraint
    forall ( i, j in TASK where i < j )
    (
        
      s[i] + d[i] <= s[j]   \/ s[j] + d[j] <= s[i]
       
    );
                  
var 0..maxt: makespan = max(e);

solve 

%:: int_search(array1d(s), input_order, indomain_random, complete)
%:: relax_and_reconstruct(array1d(b), 70)
:: restart_luby(1)
minimize makespan;
output ["makespan = ", show(makespan), "\n" ] ++
       [ show_int(3,s[t]) ++ " " 
       | t in TASK]; 
