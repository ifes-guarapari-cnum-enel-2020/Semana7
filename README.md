# Semana7
![](https://github.com/ifes-guarapari-cnum-enel-2020/Semana7/workflows/Julia%20CI/badge.svg)

Atividades Pedagógicas Não Presenciais

## Métodos de Gauss
O método de Gauss é um escalonamento de sucessivas operações elementares num sistema linear quadrático, a partir da matriz superior triangular correspondente, para resolução desse sistema linear. Os métodos de Gauss-Jacobi e Gauss-Seidel usam uma aproximação inicial da solução para esses sistemas.

https://www.ufrgs.br/reamat/CalculoNumerico/livro-py/sdsl.html

A eliminação gaussiana consiste em manipular o sistema através de determinadas operações elementares, para triangularizar a matriz estendida, e a solução é obtida por substituição regressiva. Em Julia, o método é implementado nativamente pela linguagem, em que dadas as matrizes A e B, usa-se o operador de barra invertida para achar a solução possível.
```julia
A = [1 1 1 ;
     4 4 2 ;
     2 1 -1]

B = [1 ; 2 ; 0]

X = A \ B
println(X)
```

No método de Gauss-Jacobi, cada incógnita é uma equação de ponto fixo do sistema e, a partir de aproximação inicial, é recalculada após sucessivas iterações, em que o valor anterior é usado para o próximo cálculo.
```julia
error = 10^-3

function jacobi(A, B, k)
 n = size(B,1)
 X = zeros(n)
 K = zeros(n)
 for l = 1:k
  for i = 1:n
   count = 0
   for j = 1:n
    if i != j
     count += A[i,j]*X[j]
    end
   end
   K[i] = (B[i]-count)/A[i,i]
  end
  # println(norm(X-K))
  if norm(X-K) < error
   break
  end
  X = copy(K)
 end
 return X
end

X = jacobi(A, B, 10)
println(X)
```

O método de Gauss-Seidel é uma melhoria do Jacobi, em que a próxima variável utiliza os cálculos anteriores ainda da mesma interação, diminuindo a quantidade necessária de novas repetições.
```julia
function seidel(A, B, k)
 n = size(B,1)
 X = zeros(n)
 K = zeros(n)
 for l = 1:k
  for i = 1:n
   count = [0.0 0.0]
   for j = 1:i-1
    count[1] += A[i,j]*K[j]
   end
   for j = i+1:n
    count[2] += A[i,j]*X[j]
   end
   K[i] = (B[i]-count[1]-count[2])/A[i,i]
  end
  if norm(X-K) < error
   break
  end
  X = copy(K)
 end
 return X
end

X = seidel(A, B, 100)
println(X)
```

Uma condição suficiente, porém não necessária, para que os métodos de Gauss-Seidel e Jacobi convirjam é a que a matriz seja estritamente diagonal dominante, o que deve ser observado.
