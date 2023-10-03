set i / i1*i10 /;
parameter
    f(i) / i1 1 /
    g(i) / i1 1 /
    NC /2/;
loop(i$(ord(i)>=NC),
  f(i) = f(i-2) + f(i-1);
);
display f;
