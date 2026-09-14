## 目標

證明
$$y = ReLU(W \cdot X + d)$$

$$W\in \mathbb{R}^{(3\times 3)}\quad X,b\in \mathbb{R}^{(3\times 1)}$$

$$y \in \mathbb{R}^{3\times 1}$$



## 證明分為兩個部分

$$
X,W,b
\longrightarrow
z=WX+b
\longrightarrow
y=\operatorname{ReLU}(z)
$$

1. Arithmetic Constraint System
2. Plookup ReLU



## 計算 $z = W \cdot X+b$

Define:

$$z = W \cdot X+b$$

where:

$$z = \begin{bmatrix} z_1\\ z_2\\ z_3 \end{bmatrix}$$

$$z_1 = w_{11}x_1 + w_{12}x_2 + w_{13}x_3 + b_1$$

$$z_2 = w_{21}x_1 + w_{22}x_2 + w_{23}x_3 + b_2$$

$$z_3 = w_{31}x_1 + w_{32}x_2 + w_{33}x_3 + b_3$$

## Constraint System 證明 $z = W \cdot X$

對每個輸出建立 arithmetic constraint

$$z_1 - (w_{11}x_1 + w_{12}x_2 + w_{13}x_3 + b_1) = 0$$

$$z_2 - (w_{21}x_1 + w_{22}x_2 + w_{23}x_3 + b_2) = 0$$

$$z_3 - (w_{31}x_1 + w_{32}x_2 + w_{33}x_3 + b_3) = 0$$


## 建立固定的 ReLU Lookup Table

$$z \rightarrow \operatorname{ReLU}(z)$$

定義固定的 ReLU table：

$$T_{\mathrm{ReLU}}=\{(z,y)\mid y=\operatorname{ReLU}(z)\}$$

> $\operatorname{ReLU}(z)=\max(0,z)$

Example:

| z | y |
|---:|---:|
| -3 | 0 |
| -2 | 0 |
| -1 | 0 |
| 0 | 0 |
| 1 | 1 |
| 2 | 2 |
| 3 | 3 |


## 將每個 $z_i$ 丟進 ReLU lookup

當 Constraint System 得到：
$$
Z=
\begin{bmatrix}
z_1\\
z_2\\
z_3
\end{bmatrix}
$$

定義輸出:

$$
y=
\begin{bmatrix}
y_1\\
y_2\\
y_3
\end{bmatrix}
$$

做 lookup 驗證



##  Overall implementation flow

```text
W, X, b
   │
   ▼
Arithmetic Constraint System
   │
   │ prove:
   │ Z = WX + b
   ▼
z = [z1, z2, z3]
   │
   │
   ├──── (z1, y1) ────┐
   ├──── (z2, y2) ────┼── Plookup
   └──── (z3, y3) ────┘
                       │
                       ▼
                Fixed ReLU Table
                       │
                       ▼
                y = [y1, y2, y3]