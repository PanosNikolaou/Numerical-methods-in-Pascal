{
ENVIRONMENTAL MODELLING SS 2017/2018
PANAGIOTIS NIKOLAOU
FINDING AN INVERSE MATRIX OF A 3X3 TABLE
USING THE FORMULA A^-1=1/DET(A)*ADJ(A)
SOLVE LINEAR SYSTEM BY
GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

Ax = b
roundoff error
gauss elimination
}

PROGRAM gauss;

TYPE
  TMatrix  = ARRAY[1..3, 1..3] OF Real;
  TVector  = ARRAY[1..3] OF Real;

VAR
  A,T: TMatrix;
  ADJ,AINV : TMatrix;
  B,X : TVector;
  N,i,j : Integer;
  endch: Char;
  DeterminantVAR: Real;

PROCEDURE READDATA(VAR A: TMatrix; N: Integer);

VAR
  I, J: Integer;

BEGIN
  N := 3;
  WriteLn;
  FOR I := 1 TO N DO
    BEGIN
      FOR J := 1 TO N DO
        BEGIN
          Write( I,',',J, ': ');
          ReadLn( A[I,J])
        END
    END;
  WriteLn;
  WriteLn( 'INITIAL MATRIX');
  FOR I:= 1 TO N  DO
    BEGIN
      FOR J:= 1 TO N DO
        Write( A[I,J]:7:4, '  ');
      WriteLn
    END;
  WriteLn
END;

PROCEDURE READBDATA(VAR B: TVector; N: Integer);

VAR
  I, J: Integer;

BEGIN
  N := 3;
  WriteLn;
  FOR I := 1 TO N DO
        BEGIN
          Write( I, ': ');
          ReadLn( B[I])
        END;
  //WriteLn;
  //WriteLn( 'INITIAL MATRIX');
  FOR I:= 1 TO N  DO
      BEGIN
           Write( B[I]:7:4, '  ');
           WriteLn
      END;
  WriteLn
END;

// Solve A x = b
PROCEDURE gauss (VAR A : TMatrix; VAR b,x : TVector);
VAR rowx : integer;
    i, j, k, n, m : integer;
    amax, xfac, temp, temp1 : double;
BEGIN
  n := 3;
  rowx := 0;  // Keep count of the row interchanges
  //n := A.r;

  WriteLn( 'Forward elimination started');

  FOR k := 1 TO n - 1 DO
      BEGIN
      amax := abs (A[k,k]);
      m := k;
      // Find the row with largest pivot
      FOR i := k + 1 TO n DO
          BEGIN
          xfac := abs (A[i,k]);
          IF xfac > amax THEN
             BEGIN
             amax := xfac;
             m := i;
             END;
          END;

      IF m <> k THEN
         BEGIN  // Row interchanges
         rowx := rowx+1;
         temp1 := b[k];
         b[k] := b[m];
         b[m]  := temp1;
         FOR j := k TO n DO
             BEGIN
             temp := a[k,j];
             a[k,j] := a[m,j];
             a[m,j] := temp;
             END;
      END;

      //forward elimination

      FOR i := k+1 TO n DO
          BEGIN
          xfac := a[i, k]/a[k, k];
          FOR j := k+1 TO n DO
              a[i,j] := a[i,j]-xfac*a[k,j];
          b[i] := b[i] - xfac*b[k]
          END;
      END;

 WriteLn( 'Back substitution started');

  // Back substitution
  FOR j := 1 TO n DO
      BEGIN
      k := n-j + 1;
      x[k] := b[k];
      FOR i := k+1 TO n DO
          BEGIN
          x[k] := x[k] - a[k,i]*x[i];
          END;
  x[k] := x[k]/a[k,k];
  END;
END;


PROCEDURE Cofactor (VAR TA:TMatrix;A: TMatrix);

VAR
  A11,A12,A13,A21,A22,A23,A31,A32,A33: Real;

BEGIN
  n:=3;
  A11 := A[2,2]*A[3,3] - A[2,3] * A[3,2];
  TA[1,1]:= A11;
  A12 := - (A[2,1]*A[3,3] - A[2,3] * A[3,1]);
  TA[1,2]:= A12;
  A13 := A[2,1]*A[3,2] - A[2,2] * A[3,1];
  TA[1,3]:= A13;
  A21 := - (A[1,2]*A[3,3] - A[1,3] * A[3,2]);
  TA[2,1]:= A21;
  A22 := A[1,1]*A[3,3] - A[1,3] * A[3,1];
  TA[2,2]:= A22;
  A23 := - (A[1,1]*A[3,2] - A[1,2] * A[3,1]);
  TA[2,3]:= A23;
  A31 := A[1,2]*A[2,3] - A[1,3] * A[2,2];
  TA[3,1]:= A31;
  A32 := - (A[1,1]*A[2,3] - A[1,3] * A[2,1]);
  TA[3,2]:= A32;
  A33 := A[1,1]*A[2,2] - A[1,2] * A[2,1];
  TA[3,3]:= A33;

END;

FUNCTION Determinant (A: TMatrix): Real;
VAR
  SUM: Real;

BEGIN
  SUM := A[1,1] * (A[2,2]*A[3,3] - A[3,2]*A[2,3])
       - A[1,2] * (A[2,1]*A[3,3] - A[3,1]*A[2,3])
       + A[1,3] * (A[2,1]*A[3,2] - A[3,1]*A[2,2]);
  RESULT := SUM;
END;


BEGIN

    N:=3;

    // INITIALIZE MATRIX A AND T
    BEGIN
     FOR i:= 1 TO N DO
      FOR j:= 1 TO N DO
         A[i,j]:= i * j;
         T[i,j]:= i * j;
    END;

    // INITIALIZE MATRIX B AND X
    BEGIN
     FOR i:= 1 TO N DO
         B[i]:= i;
         X[i]:= i;
    END;

    // Test matrices

    // 12  2  2 = 3
    // 3   2  3 = 2
    // 2   1  2 = 1

    //x1 = 0.1
    //x2 = 1
    //x3 = -0.1

    WriteLn( 'INSERT VALUES FOR MATRIX A');
    READDATA(A,N);

    WriteLn( 'INSERT VALUES FOR VECTOR B');
    READBDATA(B,N);

    Cofactor(T,A);

    WriteLn( 'COFACTOR MATRIX');
     FOR I:= 1 TO N  DO
       BEGIN
         FOR J:= 1 TO N DO
           Write( T[I,J]:7:4, '  ');
         WriteLn();
       END;

    DeterminantVAR := Determinant(A);
    WriteLn();
    WriteLn( 'THE DETERMINANT IS: ', DeterminantVAR:7:4);
    WriteLn();

    // IF MATRIX IS ZERO MATRIX IS NOT INVERTABLE
    IF DeterminantVAR = 0 THEN
       BEGIN
            WRITELN('MATRIX IS NOT INVERTABLE!');
            ReadLn(endch);
            EXIT();
       END;

    FOR i := LOW(A) TO HIGH(A) DO
      FOR j := LOW(A[1]) TO HIGH(A[1]) DO
        ADJ[j,i] := T[i,j];

    WriteLn( 'ADJOINT MATRIX');

    FOR I:= 1 TO N  DO
      BEGIN
        FOR J:= 1 TO N DO
          Write( ADJ[I,J]:7:4, '  ');
        WriteLn();
      END;

    FOR I:= 1 TO N  DO
        FOR J:= 1 TO N DO
          AINV[I,J] := (1/DeterminantVAR) * ADJ[I,J];
    WriteLn();

    WriteLn( 'INVERSE MATRIX');

    FOR I:= 1 TO N  DO
      BEGIN
        FOR J:= 1 TO N DO
          Write( AINV[I,J]:7:4, '  ');
        WriteLn();
      END;

    WriteLn;
    WriteLn( 'SOLVING SYSTEM AX=B BY GAUSSIAN ELIMINATION');

    gauss(a,b,x);

    WriteLn;
    WriteLn( 'MATRIX A AFTER GAUSSIAN ELIMINATION');

    FOR I:= 1 TO N  DO
        BEGIN
          FOR J:= 1 TO N DO
            Write( A[I,J]:7:4, '  ');
          WriteLn();
        END;

    WriteLn;
    WriteLn( 'X SOLUTIONS');

    FOR I:= 1 TO N  DO
        Write( X[I]:7:4, '  ');

    WriteLn;
    WriteLn;
    WriteLn( 'PRESS ENTER TO TERMINATE');
    ReadLn(endch)
END.

