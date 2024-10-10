//! exec('E:\GitHub\Facu\Scilab\fin2.sce', -1)
N = 100
A = [4 3 0; 3 4 -1; 0 -1 4];
B = [24 30 -24]'; 
// A = 8*eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-3,1),3) + diag(ones(N-3,1),-3)
// b = ones(N,1)

function s1 = remonte(A, b) // Ax = b
    n = size(A, 1)
    x = zeros(n, 1)
    x(n) = b(n)/A(n,n)
    
    for i = n-1:-1:1
        x(i) = (b(i) - A(i,i+1:n)*x(i+1:n))/A(i,i);
    end
    s1 = x
endfunction


function s1 = descenso(A, b) // Ax = b, con A triangular inferior
    n = size(A, 1)  
    x = zeros(n, 1)  
    x(1) = b(1) / A(1,1)  

    // Bucle hacia adelante
    for i = 2:n
        x(i) = (b(i) - A(i, 1:i-1) * x(1:i-1)) / A(i, i);
    end
    s1 = x  
endfunction


// Realiza la eliminación de Gauss sin pivoteo en una matriz tridiagonal. 
// Recibe una matriz A y un vector b, y devuelve la matriz triangular superior y el vector resultante.
function s2 = gauss_tridiagonal(A, b)
    // Obtener el tamaño de la matriz A
    n = size(A, 1); 
    
    // Realizar eliminación de Gauss para el resto de la matriz tridiagonal
    for k = 1:n-1
        m = A(k+1, k) / A(k, k);  
        A(k+1, :) = A(k+1, :) - (m * A(k, :));  
        b(k+1, 1) = b(k+1, 1) - (m * b(k, :));
    end
    
    s2 = b; 
    s2(n) = s2(n)/A(n,n);
    
    for i = n-1:-1: 1
        s2(i) = (b(i) - (A(i, i+1)*s2(i+1))) / A(i, i);
    end
        
endfunction



// Resuelve Ax=B por metodo de Gauss. B es una matriz. Si B=I del tamaño de A, entonces calcula la inversa.
function [A, B] = gauss_solve(A, B)
    // Obtener el tamaño de A
    n = size(A, 1);  // Suponemos que A es una matriz cuadrada de tamaño n x n

    // Aplicar eliminación de Gauss-Jordan para reducir [A | B] a [I | X]
    for k = 1:n
        pivot = A(k, k);
        A(k, :) = A(k, :) / pivot;
        B(k, :) = B(k, :) / pivot;

        // Eliminar todos los elementos por encima y por debajo del pivote en la columna k
        for i = 1:n
            if i ~= k then
                factor = A(i, k);
                A(i, :) = A(i, :) - (factor * A(k, :));
                B(i, :) = B(i, :) - (factor * B(k, :));
            end
        end
    end
    // Extraer la parte derecha de la matriz aumentada, que es la matriz transformada X
    
endfunction


function [s1, s2] = gauss(A, b)
    a = size(A)
    n = a(1)
    
    for i = 1:(n-1)
        for j = (i+1):a(1)
            mjk = A(j,i)/A(i,i)
            A(j,i)=0
            A(j,i+1:n) = A(j,i+1:n) - mjk*A(i,i+1:n)
            b(j) = b(j) - mjk*b(i)
        end
    end
    
    s1 = A
    s2 = b
endfunction


function [s1, s2] = gauss_recursivo(A, b)
    dim = size(A);  // Obtiene las dimensiones de la matriz A
    
    if dim(1) <= 1 then
        s1 = A;
        s2 = b;
        return
    end
 
    for i = 2:dim(1)
        m = (A(i, 1) / A(1, 1));  
        A(i, :) = A(i, :) - (m * A(1, :));  
        b(i, 1) = b(i, 1) - (m * b(1, 1));  
    end

    [s1_sub, s2_sub] = gauss_recursivo(A(2:$, 2:$), b(2:$));  
    
    s1 = [A(1, :); [zeros(dim(1)-1, 1), s1_sub]];
    s2 = [b(1, :); s2_sub];
endfunction


// Realiza la eliminación de Gauss sin Pivoteo. Recibe una matriz A
function [L, U] = gauss_sin_pivoteo(A)
    U = A;                      // U se inicializa como la matriz A
    m = size(A, 1);             // m es el número de filas de A
    L = eye(m,m);               // L se inicializa como la matriz identidad de tamaño mxm

    // Iteración sobre las columnas principales de la matriz para la eliminación Gaussiana
    for k = 1:m-1
        // Verifica que el pivote (U(k,k)) no sea cero para evitar divisiones por cero
        if U(k,k) == 0 then
            error("Eliminación de Gauss fallida: pivote cero encontrado en U(" + string(k) + "," + string(k) + ")");
        end

        // Realiza la eliminación de Gauss para anular los elementos debajo del pivote
        for j = k+1:m
            L(j,k) = U(j,k) / U(k,k);               // Calcula el multiplicador para la fila j
            U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m);  // Resta el múltiplo de la fila k de U a la fila j
        end
    end
endfunction


// Realiza la eliminación de Gauss con Pivoteo Parcial. Recibe A matriz
function [P,L,U] = egpp(A)
    U = A;                      // U inicializa como la matriz A
    m = size(A, 1);             // m es el número de filas de A
    L = eye(m,m);               // L es la matriz identidad de tamaño mxm
    P = eye(m,m);               // P es la matriz de permutación (identidad inicialmente)

    // Iteración sobre las columnas principales de la matriz para la eliminación Gaussiana
    for k = 1:m-1
      // Encuentra el índice del elemento máximo en valor absoluto en la columna k, para pivoteo parcial
      ind = find(abs(U(k:m,k)) == max(abs(U(k:m,k))),1); 
      ind = ind + (k-1);        // Ajusta el índice para el rango actual

      // Intercambia filas en U para llevar el pivote a la posición correcta
      U([k ind],k:m) = U([ind k],k:m); 

      // Intercambia las filas en L para mantener consistencia en la factorización
      L([k ind],1:k-1) = L([ind k],1:k-1); 

      // Intercambia las filas en P para registrar la permutación
      P([k ind],:) = P([ind k],:);

      // Realiza la eliminación de Gauss para anular los elementos debajo del pivote
      for j = k+1:m
          L(j,k) = U(j,k) / U(k,k);                 // Calcula el multiplicador para la fila j
          U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m);    // Resta el múltiplo de la fila k de U a la fila j
      end 
    end
endfunction


// Resuelve un sistema de ecuaciones lineales Ax = b con Gauss y pivoteo parcial.
function x = Gesolver(A, b)
    [P, L, U] = egpp(A);

    // Obtiene el tamaño de la matriz `L`, que corresponde al número de ecuaciones
    sz = size(L, 1);

    // Ajusta el vector `b` de acuerdo con la matriz de permutación `P`.
    b = P * b;

    // Inicializa el vector `c` para la sustitución hacia adelante.
    c(1) = b(1);

    // Sustitución hacia adelante para resolver el sistema Lc = Pb.
    for i = 2:sz
        suma = 0;
        // Calcula la suma de L(i, j) * c(j) para las componentes anteriores.
        for j = 1:i-1
            suma = suma + c(j) * L(i, j);
        end
        // Calcula el valor de c(i) usando la fórmula: c(i) = b(i) - suma
        c(i) = b(i) - suma;
    end

    // Inicializa la última componente de `x` para la sustitución hacia atrás.
    x(sz) = c(sz) / U(sz, sz);

    // Sustitución hacia atrás para resolver el sistema Ux = c.
    for i = 1:sz-1
        suma = 0;
        // Calcula la suma de U(sz-i, sz-j+1) * x(sz-j+1) para las componentes anteriores.
        for j = 1:i
            suma = suma + x(sz-j+1) * U(sz-i, sz-j+1);
        end
        // Calcula el valor de x(sz-i) usando la fórmula: x(i) = (c(i) - suma) / U(i, i)
        x(sz-i) = (c(sz-i) - suma) / U(sz-i, sz-i);
    end
endfunction


//Método factorización Doolittle. Recibe matriz A 
function [L,U] = dolittle(A)
    sz = size(A,1)
    L = eye(sz,sz)
    U = zeros(sz,sz)
    
    U(1,:) = A(1,:)
    for k = 1:sz
           
        for j = k:sz
            suma = 0
            for m = 1:k-1
                suma = suma + L(k,m)*U(m,j)
            end
            U(k,j)=A(k,j) - suma
            
            for i = k+1:sz
                suma = 0
                for m = 1:k-1
                    suma = suma + L(i,m)*U(m,k)
                end    
                L(i,k)=(A(i,k) - suma)/U(k,k)
            end
        end
    end
endfunction

//Resuelve un sistema con factorización Doolittle. Recibe matriz A y vector b
function x = dlsolver(A,b)
    sz = size(A,1)
    [L,U] = dolittle(A)
    c(1) = b(1)
    
    for i = 2:sz
        suma = 0
        for j = 1:i-1
           suma = suma + c(j)*L(i,j)   
        end
        c(i) = b(i)-suma 
    end
    
    x(sz) = c(sz)/U(sz,sz)

    for i = 1:sz-1
        suma = 0
        for j = 1:i
           suma = suma + x(sz-j+1)*U(sz-i,sz-j+1)   
        end
        x(sz-i) = (c(sz-i)-suma)/U(sz-i,sz-i) 
    end    
    
endfunction

//Método de factorización Cholesky. Recibe matriz A.
function [U, ind] = Cholesky(A)
    eps = 1.0e-8
    n = size(A,1)
    U = zeros(n,n)
    iter = 0;

    for k = 1:n
        t = A(k,k) - U(1:k,k)'*(U(1:k,k))
        if (t <= eps)
            mprintf("Matriz no definida positiva.\n")
            ind = 0
            return
        end
        if (A <> A')
            mprintf("Matriz no simetrica. \n")
            ind = 0
            return 
        end
        U(k,k)= sqrt(t)
        for j = k+1:n
            U(k,j) = ( A(k,j) - U(1:k,k)'*U(1:k,j) )/U(k,k)
            iter = iter+ 1;
        end
    end
    ind = 1
    mprintf("Cantidad de iteraciones: %d\n", iter);
endfunction

//Resuelve un sistema con factorización Cholesky. Recibe matriz A y vector b
function x = chsolver(A,b)
    [U,ind] = Cholesky(A)
    if (ind == 1)
        c(1) = b(1)
        L = U'
        
        for i = 2:sz
            suma = 0
            for j = 1:i-1
               suma = suma + c(j)*L(i,j)   
            end
            c(i) = b(i)-suma 
        end
        
        x(sz) = c(sz)/U(sz,sz)
    
        for i = 1:sz-1
            suma = 0
            for j = 1:i
               suma = suma + x(sz-j+1)*U(sz-i,sz-j+1)   
            end
            x(sz-i) = (c(sz-i)-suma)/U(sz-i,sz-i) 
        end  
    else
        mprintf("Error, no se puede \n");
    end
endfunction

// ===============================================================================
// ===============================================================================
// ===============================================================================


function [A_pivoted, row_exchange] = pivoteo_parcial(A, k)
    // A: matriz de entrada
    // k: columna en la que se realiza el pivoteo (1-indexado en Scilab)
    // A_pivoted: matriz con las filas intercambiadas tras el pivoteo
    // row_exchange: matriz de intercambio de filas, mostrando qué filas se intercambiaron

    [n, m] = size(A);
    
    // Verificamos que el índice de pivote esté dentro de los límites de la matriz
    if k > n then
        error("El índice del pivote está fuera del rango");
    end

    row_exchange = eye(n, n);
    
    // Encontrar la fila con el mayor valor absoluto en la columna k desde la fila k en adelante
    [max_val, max_row] = max(abs(A(k:n, k)));
    
    // max_row nos da el índice relativo (empezando desde k), así que ajustamos
    max_row = max_row + k - 1;
    
    // Si la fila con el mayor pivote no es la actual, intercambiamos las filas
    if abs(A(max_row, k)) > abs(A(k, k)) then
        // Intercambiamos las filas en la matriz
        temp_row = A(k, :);          
        A(k, :) = A(max_row, :);     
        A(max_row, :) = temp_row;   
        
        // También intercambiamos las filas en la matriz de intercambio
        temp_exchange = row_exchange(k, :);
        row_exchange(k, :) = row_exchange(max_row, :);
        row_exchange(max_row, :) = temp_exchange;
    end

    A_pivoted = A;
endfunction


//Test de Jacobi para saber si converge a la solución. Recibe matriz A
function radspec = jacobitest(A)
    Nor = NOR_jacobi(A)
    radspec = norm(eigs(Nor), 'inf')
    
    if(radspec < 1)
        mprintf("Converge para cualquier x0 inicial.\n")
    else
        mprintf("No converge para cualquier x0 inicial.\n")
    end
endfunction


//Método de resolución de Jacobi. Recibe matriz A, vector b, un x0 inicial y un épsilon.
function xn = jacobisolver(A,b,x0,eps) 
    n = size(A,1)
    xn = x0
    xk = x0
    suma = 0
    cont = 0
    
    while(norm(xn - xk) > eps || cont == 0) 
        xk = xn;
        for i = 1:n
            suma = 0;
            for j = 1:n 
                if (i <> j)
                    suma = suma + A(i,j)*xk(j)
                end
            end
            xn(i) = 1/(A(i,i))*(b(i)-suma)
        end
     cont = cont + 1;
    end
    
    mprintf("Cantidad de iteraciones: %d\n",cont);
endfunction


//Test de Gauss-Seidel para saber si converge a la solución. Recibe matriz A
function radspec = gausstest(A)
    Nor = NOR_gaussSeidel(A);    
    radspec = norm(eigs(Nor), 'inf')
    
    if(radspec < 1)
        mprintf("Converge para cualquier x0 inicial.\n")
    else
        mprintf("No converge para cualquier x0 inicial.\n")
    end
endfunction


//Método de resolución de Gauss-Seidel. Recibe matriz A, vector b, un x0 inicial y un épsilon.
function xn = gausssolver(A,b,x0,eps)  
    n = size(A,1)
    xn = x0
    xk = x0
    suma = 0
    cont = 0
    
    while(abs(norm(xn-xk)) > eps | cont == 0) 
        xk = xn
        for i = 1:n
            suma = 0
            for j = 1:i-1 
                suma = suma + A(i,j)*xn(j)
            end
            
            for j = i+1:n
                suma = suma + A(i,j)*xn(j)
            end
            xn(i) = 1/(A(i,i))*(b(i)-suma)
        end
     cont = cont + 1
    end
    
    mprintf("Cantidad de iteraciones: %d\n",cont);
endfunction


//Factorización Gram Schmidt. Recibe matriz A
function GS = Gram_Schmidt(A)   // A debe tener columnas LI
    sz = size(A,2)
    Qu(:,1) = A(:,1)/norm(A(:,1))
    for i = 2:sz
        suma = 0
        for j = 1:i-1
            suma = suma + (A(:,j)'*Qu(:,j))*Qu(:,j)
        end
        Qu(:,i) = A(:,i) - suma
        Qu(:,i) = Qu(:,i)/norm(Qu(:,i))  
    end
endfunction

//Factorización QR. Recibe una matriz A.
function [Q,R] = FactoRQ(A)   // A debe tener columnas LI
    sz = size(A,2)
    Q(:,1) = A(:,1)/norm(A(:,1))
    V(1) = norm(A(:,1))
    for i = 2:sz
        suma = 0
        for j = 1:i-1
            suma = suma + (A(:,i)'*Q(:,j))*Q(:,j)
        end
        Q(:,i) = A(:,i) - suma
        V(i) = norm(Q(:,i))
        Q(:,i) = Q(:,i)/V(i)  
    end
    
    R = diag(V)
    
    
    for i = 1:sz
        for j = i+1:sz
           R(i,j) = A(:,j)'*Q(:,i) 
        end
    end
endfunction


function k = condicionMatriz(A)
    if det(A) ~= 0
        k = norm(A)* norm(inv(A));
    else 
        error("El determinante es 0. Pivotear A antes de aplicar la funcion con pivoteo_parcial.\n");
    end
endfunction


function b = is_def_pos(A)
    minEigenvalue = min(spec(A))
    if minEigenvalue > 0 then
        b = 1
    else 
        b = 0
    end
endfunction


function b = is_simetric(A)
    if A <> A' then 
        b = 0;
    else 
        b = 1
    end
endfunction


function omega = SOR_optimo(A)
    radspec = norm(eigs(NOR_jacobi(A)), 'inf');
    omega = 2 / (1 + sqrt(1 - radspec^2));
endfunction


function x = resolver(A, b)

    x = zeros(1, size(A, 1));
    [U, ind] = Cholesky(A)

    if ind == 1
        mprintf("Cholesky: \n")
        disp(U)
        g = gausssolver_L(U', b, [0 0 0], 1e-6) // R'g=b
        x = gausssolver_U(U, g, [0 0 0], 1e-6)'  // Rx = g
    else
        error("La matriz no cumple las condiciones para factorizarse por Cholesky.\n") 
    end


endfunction


function xn = jacobisolver_L(A, b, x0, eps)
    n = size(A,1)  
    xn = x0        
    xk = x0        
    suma = 0
    cont = 0       
    
    while(norm(xn - xk) > eps || cont == 0) 
        xk = xn;  // Actualizar la solución anterior
        for i = 1:n
            suma = 0;
            // Solo iterar sobre las columnas por debajo de la diagonal (triangular inferior)
            for j = 1:i-1
                suma = suma + A(i,j)*xk(j);
            end
            // Aplicar la fórmula de Jacobi para la matriz triangular inferior
            xn(i) = (b(i) - suma) / A(i,i);
        end
        cont = cont + 1;  // Incrementar el contador de iteraciones
    end
    
    // Mostrar cantidad de iteraciones
    mprintf("Cantidad de iteraciones: %d\n", cont);
endfunction


function xn = jacobisolver_U(A, b, x0, eps)
    n = size(A,1)  
    xn = x0        
    xk = x0       
    suma = 0
    cont = 0       
    
    while(norm(xn - xk) > eps || cont == 0) 
        xk = xn;  
        for i = 1:n
            suma = 0;
            // Solo iterar sobre las columnas por encima de la diagonal (triangular superior)
            for j = i+1:n
                suma = suma + A(i,j)*xk(j);
            end
            // Aplicar la fórmula de Jacobi para la matriz triangular superior
            xn(i) = (b(i) - suma) / A(i,i);
        end
        cont = cont + 1;  // Incrementar el contador de iteraciones
    end
    
    // Mostrar cantidad de iteraciones
    mprintf("Cantidad de iteraciones: %d\n", cont);
endfunction


function xn = gausssolver_L(A, b, x0, eps)  
    n = size(A, 1);    
    xn = x0;           
    xk = x0;           
    suma = 0;          
    cont = 0;          
    
    // Iteración hasta que la diferencia entre xn y xk sea menor que eps
    while (abs(norm(xn - xk)) > eps || cont == 0) 
        xk = xn;  // Actualizar la solución anterior
        for i = 1:n
            suma = 0;
            // Sumar los elementos por encima de la diagonal
            for j = 1:i-1 
                suma = suma + A(i, j) * xn(j);
            end
            
            // Sumar los elementos por debajo de la diagonal
            for j = i+1:n
                suma = suma + A(i, j) * xk(j);
            end
            
            // Actualizar la solución actual
            xn(i) = (b(i) - suma) / A(i, i);
        end
        cont = cont + 1;  // Incrementar el contador de iteraciones
    end
    
    mprintf("Cantidad de iteraciones: %d\n", cont);
endfunction


function xn = gausssolver_U(A, b, x0, eps)  
    n = size(A, 1);    
    xn = x0;           
    xk = x0;           
    suma = 0;          
    cont = 0;          
    
    // Iteración hasta que la diferencia entre xn y xk sea menor que eps
    while (abs(norm(xn - xk)) > eps || cont == 0) 
        xk = xn;  // Actualizar la solución anterior
        for i = 1:n
            suma = 0;
            // Sumar los elementos por debajo de la diagonal
            for j = 1:i-1 
                suma = suma + A(i, j) * xn(j);
            end
            
            // Sumar los elementos por encima de la diagonal
            for j = i+1:n
                suma = suma + A(i, j) * xk(j);
            end
            
            // Actualizar la solución actual
            xn(i) = (b(i) - suma) / A(i, i);
        end
        cont = cont + 1;  // Incrementar el contador de iteraciones
    end
    
    mprintf("Cantidad de iteraciones: %d\n", cont);
endfunction


function N = NOR_jacobi(A)
    sz = size(A,1)
    N = eye(sz,sz)
    I = N

    for i = 1:sz

        if A(i, i) == 0 
            [A, P] = pivoteo_parcial(A, i);
        end

        N(i,i) = N(i,i) * A(i, i)
    end

    N = (I-inv(N)*A)    
endfunction


function N = NOR_gaussSeidel(A)
    [n, m] = size(A)
    I = eye(n,n)
    N = zeros(n,n)
     
    for i = 1:n
        for j = 1:m

            if A(i, i) == 0 
                [A, P] = pivoteo_parcial(A, i);
            end

            if(i >= j)
                N(i,j) = A(i,j)
            end
        end
    end
    
    N = (I-inv(N)*A)    

endfunction


function xn = gauss_SOR(A, b, x0, eps, SOR)  
    // Definir el tamaño de la matriz A
    n = size(A, 1);
    xn = x0;
    xk = x0;
    suma = 0;
    cont = 0;

    // Iterar hasta que la norma de la diferencia sea menor que el umbral 'eps'
    while (abs(norm(xn - xk)) > eps | cont == 0) 
        xk = xn; // Guardar la iteración anterior

        for i = 1:n
            suma = 0;

            // Sumar los términos anteriores
            for j = 1:i-1 
                suma = suma + A(i,j) * xn(j);
            end

            // Sumar los términos posteriores
            for j = i+1:n
                suma = suma + A(i,j) * xn(j);
            end

            // Actualizar con sobrerelajación
            xn(i) = (1 - SOR) * xn(i) + SOR * (1 / A(i,i)) * (b(i) - suma);
        end

        cont = cont + 1; // Incrementar el contador de iteraciones
    end

    // Imprimir la cantidad de iteraciones realizadas
    mprintf("Cantidad de iteraciones: %d\n", cont);
endfunction

