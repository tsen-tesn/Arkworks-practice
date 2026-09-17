# Arkworks 練習：用 Plookup 驗證 Integer ReLU

本練習使用 Arkworks 0.5 實作 ReLU lookup proof。測試案例分為整數與浮點數兩類；現階段只實作 Integer ReLU，浮點數案例保留為後續擴充。目標不是在有限域中直接比較大小，而是證明每一組 witness row `(x, y)` 都存在於公開的 ReLU table：

$$
y=\text{ReLU}(x)=\max(0,x).
$$

實作流程參考 [kevaundray/plookup](https://github.com/kevaundray/plookup) 的 4-bit XOR 範例。該範例將 `(left, right, xor)` 三欄壓縮後做 lookup；本練習改為將 ReLU 的 `(input, output)` 兩欄壓縮，其餘 sorted concatenation、grand product、quotient polynomial 與 KZG opening 的概念相同。

> 這是一個用來理解論文協定的 POC，不是完整的 PLONK proving system。第一階段先完成並直接檢查 Plookup identities，第二階段再加入 polynomial commitment。

## 1. 建立公開 ReLU table（`src/lookup/table/relu.rs`、`src/lookup/table/mod.rs`）

將輸入範圍限制為：

$$
x\in[-32,31].
$$

公開 table 共有 64 列：

$$
T=\{(x,\max(0,x))\mid -32\le x\le31\}.
$$

```text
(-32, 0)
(-31, 0)
...
(-1, 0)
(0, 0)
(1, 1)
...
(31, 31)
```

本練習使用 BLS12-381 scalar field `Fr`。有號整數的編碼方式為：

```text
 15 -> Fr::from(15u64)
  0 -> Fr::from(0u64)
-24 -> -Fr::from(24u64)
```

有限域本身沒有正負號或大小順序；`-24` 只是 `24` 的加法反元素。因此，ReLU 的語意完全由公開 table 的合法 rows 定義，不能在 field element 上直接呼叫 `max()` 或做一般整數排序。

## 2. 建立 witness rows（`src/lookup/lookup.rs`）

題目包含兩組輸入序列。

### 整數案例（目前實作）

輸入為：

```rust
let inputs: [i64; 9] = [2, 4, -1, 0, 15, 6, -24, -2, 15];
```

對應的輸出為：

```rust
let expected: [i64; 9] = [2, 4, 0, 0, 15, 6, 0, 0, 15];
```

所以有效 witness rows 是：

```text
(2, 2), (4, 4), (-1, 0), (0, 0), (15, 15),
(6, 6), (-24, 0), (-2, 0), (15, 15)
```

### 浮點數案例（後續擴充）

輸入為：

```rust
let inputs = [2.0, 4.0, -1.0, 0.0, 0.85, 0.3, 24.0, -0.2, 0.15];
```

預期輸出為：

```rust
let expected = [2.0, 4.0, 0.0, 0.0, 0.85, 0.3, 24.0, 0.0, 0.15];
```

這組案例現階段不實作，也不納入 proof 測試。有限域不能直接表示 `f32` 或 `f64`；後續階段必須先決定 fixed-point scale、可表示範圍與量化規則，再建立相應的 ReLU lookup table。例如選定 scale $S$ 後，以整數 $\widetilde{x}=\mathrm{round}(Sx)$ 表示浮點數，並在 proof 外部將結果解碼回 $\widetilde{x}/S$。

令 evaluation domain 的大小為 $N=64$，並令 $n=N-1=63$。論文的基本形式使用長度 $N$ 的 table 與長度 $n$ 的 witness，因此用合法 row `(0, 0)` 將 witness 補到 63 列。Padding 也是 proof statement 的一部分，所以只能使用 table 中存在的 row。

## 3. 用隨機挑戰壓縮每一列（`src/lookup/proof.rs`、`src/multiset/multiset.rs`）

Plookup 的核心包含關係處理的是單一 field element。Verifier 從 transcript 取得隨機挑戰 $\theta$，將 ReLU 的兩欄綁成一個值：

$$
t_i=T_{X,i}+\theta T_{Y,i},
$$

$$
f_i=X_i+\theta Y_i.
$$

接下來要證明的是多重集合包含關係：

$$
f\subseteq t.
$$

必須以同一個 $\theta$ 同時壓縮 table 與 witness。這樣 prover 不能分別挑一個合法 input 和另一列的合法 output，再把兩者拼成非法的 ReLU row。除了隨機壓縮發生碰撞的 negligible probability 外，非法 row 不會等於任何合法 table row。

> 參考 repo 的 XOR table 有三欄，因此使用形如 $a+\theta b+\theta^2c$ 的壓縮；ReLU 只有兩欄，所以使用 $x+\theta y$。

## 4. 依 table 順序建立 sorted concatenation（`src/multiset/multiset.rs`）

將壓縮後的 witness 與 table 合併，再依「公開 table 的 row 順序」排列：

$$
s=\mathrm{sort}_t(f\mathbin\Vert t).
$$

這不是對 field representatives 做數值排序。實作時應保留每個 witness row 對應的 table index，或建立 `row -> table_index` 的 mapping，然後按 index 排序。

若某個 table row 在 witness 中出現 $k$ 次，它在 $s$ 中會連續出現 $k+1$ 次：一次來自 table，另外 $k$ 次來自 witness。這使 $s$ 的相鄰 pair 只有兩種形式：

- `(v, v)`：由 witness 的重複值產生；
- `(t_i, t_{i+1})`：沿著公開 table 前進。

長度關係為：

```text
len(f) = n      = 63
len(t) = n + 1  = 64
len(s) = 2n + 1 = 127
```

如果 witness row 不存在於 table，建立 sorted concatenation 時就應回傳錯誤；測試不可只停在這個 host-language 檢查，還要直接注入錯誤的 `f`，確認後續 grand-product identity 也會失敗。

## 5. 將 vectors 插值成 polynomials（`src/multiset/multiset.rs`、`src/multiset/multiset_equality.rs`）

建立大小 $N=64$ 的 radix-2 multiplicative subgroup：

$$
H=\{1,\omega,\omega^2,\ldots,\omega^n\},\qquad \omega^N=1.
$$

用 inverse FFT／Lagrange interpolation 建立：

- $f(X)$：63 個 witness values，最後一個 domain slot 不參與 transition；
- $t(X)$：64 個 table values；
- $h_1(X)$、$h_2(X)$：共同承載 127 個 sorted values；
- $Z(X)$：grand-product accumulator。

將 $s=(s_0,\ldots,s_{2n})$ 拆成兩組各 64 個 evaluations：

$$
h_1(\omega^i)=s_i,\qquad 0\le i\le n,
$$

$$
h_2(\omega^i)=s_{n+i},\qquad 0\le i\le n.
$$

兩段刻意重疊一個元素，因此必須滿足：

$$
h_1(\omega^n)=h_2(1).
$$

## 6. 建立 grand-product accumulator（`src/multiset/multiset_equality.rs`）

Verifier 再從 transcript 取得隨機挑戰 $\beta,\gamma\in\mathbb F$。對每個 $0\le i<n$，定義：

$$
\mathrm{NUM}_i=(1+\beta)(\gamma+f_i)
\bigl(\gamma(1+\beta)+t_i+\beta t_{i+1}\bigr),
$$

$$
\mathrm{DEN}_i=
\bigl(\gamma(1+\beta)+h_{1,i}+\beta h_{1,i+1}\bigr)
\bigl(\gamma(1+\beta)+h_{2,i}+\beta h_{2,i+1}\bigr).
$$

從 $Z(1)=1$ 開始逐列累積：

$$
Z(\omega^{i+1})=Z(\omega^i)\frac{\mathrm{NUM}_i}{\mathrm{DEN}_i}.
$$

$\beta$ 將相鄰 values 壓縮成一個 ordered pair；$\gamma$ 將 multiset equality 隨機化。若 $f\subseteq t$，sorted sequence 中的 `(v,v)` 與 `(t_i,t_{i+1})` 恰好會重組成 numerator 中相同的 factors，因此：

$$
Z(\omega^n)=\prod_{i=0}^{n-1}\frac{\mathrm{NUM}_i}{\mathrm{DEN}_i}=1.
$$

實作 inversion 前要檢查 denominator 是否為零；toy protocol 可以重新抽 challenge，Fiat–Shamir 版本則應把失敗視為 proof generation failure，不可默默用零代替 inverse。

## 7. 檢查四組 polynomial constraints（`src/multiset/quotient_poly.rs`）

令 $L_0(X)$ 與 $L_n(X)$ 分別為 domain 第一點與最後一點的 Lagrange basis polynomial。需要檢查：

### 起點

$$
L_0(X)(Z(X)-1)=0.
$$

這保證 $Z(1)=1$。

### Grand-product transition

$$
\begin{aligned}
&(X-\omega^n)Z(X)(1+\beta)(\gamma+f(X))
\bigl(\gamma(1+\beta)+t(X)+\beta t(\omega X)\bigr)\\
={}&(X-\omega^n)Z(\omega X)
\bigl(\gamma(1+\beta)+h_1(X)+\beta h_1(\omega X)\bigr)\\
&\qquad\cdot
\bigl(\gamma(1+\beta)+h_2(X)+\beta h_2(\omega X)\bigr).
\end{aligned}
$$

因子 $(X-\omega^n)$ 關閉最後一列，避免 transition 從 domain 終點繞回起點。

### $h_1/h_2$ 銜接

$$
L_n(X)\bigl(h_1(X)-h_2(\omega X)\bigr)=0.
$$

在 $X=\omega^n$ 時，$\omega X=1$，所以這正是 $h_1(\omega^n)=h_2(1)$。

### 終點

$$
L_n(X)(Z(X)-1)=0.
$$

這保證完整乘積回到 1，也就是 numerator 與 denominator 所描述的 multisets 相等。

Phase A 應先在 $H$ 的每一點直接檢查上述 constraints。這一步最容易發現 off-by-one、rotation 與 $h_1/h_2$ overlap 錯誤。

## 8. 建立 constraint 與 quotient polynomial（`src/multiset/quotient_poly.rs`、`src/multiset/proof.rs`）

進入 committed protocol 後，使用新的 transcript challenge $\eta$ 隨機合併四組 constraints，避免不同 constraints 的錯誤互相抵消。令合併後的 polynomial 為 $C(X)$。

Domain 的 vanishing polynomial 是：

$$
Z_H(X)=X^N-1=X^{64}-1.
$$

若所有 constraints 都在 $H$ 上成立，則 $C(X)$ 可被 $Z_H(X)$ 整除：

$$
Q(X)=\frac{C(X)}{Z_H(X)},
$$

且 polynomial division 的 remainder 必須為零。實作時不要只計算某一點上的 $C(\zeta)/Z_H(\zeta)$；prover 必須先由完整的 $C(X)$ 建立並承諾真正的 quotient polynomial $Q(X)$。

## 9. KZG commitments 與 openings（`src/kzg10.rs`、`src/transcript.rs`、`src/multiset/proof.rs`）

Phase B 使用 KZG 將 polynomial identities 轉成 succinct proof：

1. 對公開 table 的兩個 columns 做 preprocess／commitment。
2. Prover 承諾 witness 的 input 與 output columns；transcript 必須先有不可由 prover 事後調整的資料。
3. Transcript 產生 $\theta$，prover 與 verifier 利用 commitment 的線性性得到壓縮後的 $f(X)$、$t(X)$ 及其 commitments。
4. Prover 建立 sorted concatenation，並承諾 $h_1(X)$、$h_2(X)$。
5. Transcript 吸收上述 commitments 後產生 $\beta,\gamma$。
6. Prover 建立並承諾 $Z(X)$。
7. Transcript 產生 constraint-combination challenge $\eta$。
8. Prover 建立並承諾 $Q(X)$。
9. Transcript 產生 evaluation challenge $\zeta\notin H$。
10. Prover 提供 verifier 計算 constraints 所需的 evaluations：

$$
\begin{gathered}
f(\zeta),\ t(\zeta),\ t(\omega\zeta),\\
h_1(\zeta),\ h_1(\omega\zeta),\\
h_2(\zeta),\ h_2(\omega\zeta),\\
Z(\zeta),\ Z(\omega\zeta),\ Q(\zeta).
\end{gathered}
$$

11. Verifier 先用 KZG opening proofs 確認 evaluations 來自已承諾的 polynomials，再檢查：

$$
C(\zeta)=Z_H(\zeta)Q(\zeta).
$$

同一 evaluation point 的 openings 可以 batch；$\zeta$ 與 $\omega\zeta$ 是兩個不同的 opening points。Fiat–Shamir transcript 必須吸收 protocol label、公開 table／domain、所有 commitments 與先前訊息，且 challenge 的產生順序不可交換。參考 POC 特別提醒：若在 transcript 尚為空時就抽取 row-compression challenge，non-interactive 版本可能缺乏足夠 entropy；因此本練習應先承諾原始 witness columns，或把 lookup 嵌入已有 transcript state 的外層 protocol。

## 10. 建議的實作順序（`src/lookup/`、`src/multiset/`、`src/kzg10.rs`、`src/transcript.rs`）

先完成可直接測試的 Phase A，再加入 Phase B：

```text
ReLU rows 與 signed-field encoding
              ↓
witness 建立與合法 padding
              ↓
用 θ 壓縮 table / witness rows
              ↓
依 table index 建立 sorted concatenation
              ↓
插值 f、t、h₁、h₂
              ↓
用 β、γ 建立 Z
              ↓
逐點檢查四組 constraints
              ↓
建立 C 與 quotient Q，確認 remainder = 0
              ↓
加入 transcript、KZG commitments 與 openings
```

### 與參考 repo 對齊的檔案結構

以下沿用 `kevaundray/plookup` 的實際目錄分層，只把 XOR 專用的 `four_bits.rs` 換成 ReLU 專用的 `relu.rs`：

```text
plookup/
├── Cargo.toml
├── README.md
├── src/
│   ├── lib.rs
│   ├── kzg10.rs
│   ├── transcript.rs
│   ├── lookup/
│   │   ├── mod.rs
│   │   ├── lookup.rs
│   │   ├── proof.rs
│   │   └── table/
│   │       ├── mod.rs
│   │       ├── generic.rs
│   │       └── relu.rs
│   └── multiset/
│       ├── mod.rs
│       ├── multiset.rs
│       ├── multiset_equality.rs
│       ├── proof.rs
│       └── quotient_poly.rs
└── tests/
    └── lookup.rs
```

## 11. 最低測試集合（`tests/lookup.rs`）

### 現階段：整數案例

- 公開 table 恰好包含 `x = -32..31` 的 64 個 ReLU rows。
- 整數輸入 `[2, 4, -1, 0, 15, 6, -24, -2, 15]` 得到 `[2, 4, 0, 0, 15, 6, 0, 0, 15]`，padding 後 witness 長度為 63。
- 重複查詢（例如 `(15, 15)`）在 sorted concatenation 中保留正確 multiplicity。
- $h_1$ 與 $h_2$ 各有 64 個 evaluations，且交界值相同。
- 合法 witness 的 $Z$ 起點與終點都是 1，四組 constraints 在整個 $H$ 上為零。
- 非法輸出（例如 `(2, 3)`）、超出範圍的輸入（例如 `(32, 32)`）與被竄改的 compressed witness 都無法通過。
- 合法 case 的 quotient division remainder 為零；非法 case 的 constraints 或 remainder 非零。
- Phase B 中竄改 evaluation、opening proof、commitment 或 transcript message 都必須驗證失敗。

### 後續階段：浮點數案例

- 輸入 `[2.0, 4.0, -1.0, 0.0, 0.85, 0.3, 24.0, -0.2, 0.15]` 得到 `[2.0, 4.0, 0.0, 0.0, 0.85, 0.3, 24.0, 0.0, 0.15]`。
- 測試 fixed-point encode/decode、量化誤差、邊界值與超出可表示範圍的輸入。
- 浮點數測試在 fixed-point 規格確定前標記為尚未實作，不使用原生浮點數直接建立 field witness。

## 參考資料

- [plookup: A simplified polynomial protocol for lookup tables](https://eprint.iacr.org/2020/315)
- [kevaundray/plookup：4-bit XOR POC](https://github.com/kevaundray/plookup)
