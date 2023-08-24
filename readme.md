# 整数质因数分解 N = 21 —— 基于 *Shor* 算法的实现

## 1. 基本介绍

本题的基本问题是针对 `N = 21` 的情况，实现基于量子算法的质因数分解。正如题目中也提示到的，Shor 算法正是一种量子质因数分解算法，其相对于经典算法具有指数级的理论加速，应用前景被极大的看好，甚至认为威胁到了经典通讯领域的 RSA 加密算法。本组选择基于 pyqpanda 实现 Shor 算法以求解 `N = 21` 的质因数分解。


### 1.1 关于 Shor 算法

Shor 算法的相关内容已有很多的资料，此处不详细展开讨论。仅给一个简单的介绍：

Shor 算法主要思想在于将分解质因数问题转化为寻找模指数电路的周期问题，相关的步骤利用到了一些数论中的定理。直观的来说，以 `N = 21` 与 `a = 2` 为例，`f(x) = a ^ x mod N` 将会是一个周期函数：

```
f(1) = 2 ^ 1 mod 21 = 2
f(2) = 2 ^ 2 mod 21 = 4
f(3) = 2 ^ 3 mod 21 = 8
f(4) = 2 ^ 4 mod 21 = 16
f(5) = 2 ^ 5 mod 21 = 11
f(6) = 2 ^ 6 mod 21 = 1
f(7) = 2 ^ 7 mod 21 = 2
f(8) = 2 ^ 8 mod 21 = 4
...
```

对上面的例子，可求得周期 `r = 6`，则 `a ^ (r / 2) = 8`，Shor 算法利用数论的知识表明，此时 `a ^ (r / 2) + 1` 与 `a ^ (r / 2) - 1` 关于 `N = 21` 的最大公因数，即两个对应的 `N = 21` 的质因数。具体来说，对于上述的例子有 `gcd(a ^ (r / 2) + 1, N) = gcd(9, 21) = 3` 与 `gcd(a ^ (r / 2) - 1, N) = gcd(7, 21) = 7`。从而得到了 `N = 21` 的一对质因数。

可见，如何寻找 `f(x) = a ^ x mod N` 的周期将是 Shor 算法要解决的重点。Shor 算法使用一个常称为 `period-finding` 的量子子程序来求解。而这个部分的量子子程序往往需要构造相应的量子线路，其中的模指数电路（Fig. 1 中间部分）的有效设计与实现被认为是 Shor 算法的一大瓶颈。

<br />
<div style="margin: 0 auto; width: 760px;
text-align: center;">
    <img src="./images/1. period finding.jpg" style="margin-bottom:17px"/>
    <label style="font-size: 17px;"><b>Fig. 1</b> period finding 量子线路基本构造</label>
</div><br />


需要一提的是，题目限定了 `N = 21`，但对于所选择的 `a` 并未限定。事实上如果针对特定的 `a` 并辅以先验信息，可以设计出十分简洁（特化）的 period finding 线路。例如 [<sup>1</sup>]提出的如下实现（Fig. 2）做到了只需用 5 个量子比特就能构建并求解。但相应的是，这样的设计运用到了大量的先验信息，对于这里的例子而言，是考虑到了：1. 已知针对 `a=4, N=21` 的例子，有周期 `r=3`。2. 且对于 `a=4, N=21` 计算可知，其 `f(x) = a ^ x mod N` 的函数值仅取 1、4、16 三个 2 的 n 次幂值（事实上，原文给出的优化思路将之进一步看作 4 的 n 次幂值，从而可以用量子比特 |00⟩ 表示 4 的 0 次方 1，用量子比特 |01⟩ 表示 4 的 1 次方 4，用量子比特 |10⟩ 表示 4 的 2 次方 16，从而仅需两个量子比特就可以表征 `f(x) = a ^ x mod N` 的运算结果）。

<br />
<div style="margin: 0 auto; width: 700px;
text-align: center;">
    <img src="./images/2. a=4.png" style="margin-bottom:17px"/>
    <label style="font-size: 17px;"><b>Fig. 2</b> 先验辅助设计的 a=4, N=21 period finding 量子线路</label>
</div><br />

但显然过分依赖先验来完成 Shor 算法在通用性上有较大的缺陷，因此本组没有采用上面的方案，实现目标是针对例子 `N = 21` 的所有合法的 `a` 的取值均可以展开计算，最终我们将使用 8 个量子比特的实现来求解其中的周期。这部分的详细讨论见后文。

### 1.2 关于项目

为了便于理解整个的流程，此处给出项目的基本结构，见下图 Fig. 3。实现基于 python(3.7.0)。
<br />
<div style="margin: 0 auto; width: 990px;
text-align: center;">
    <img src="./images/3. process.png" style="margin-bottom:17px"/>
    <label style="font-size: 17px;"><b>Fig. 3</b> 项目结构与运行流程</label>
</div><br />

<br />


## 2. 程序实现

### 2.1 主要流程

该 Shor 算法的实现参照 <a href="http://mmrc.amss.cas.cn/tlb/201702/W020170224608149940643.pdf"><cite>Quantum Computation and Quantum Information</cite>(p 233-234)</a>[<sup>2</sup>] 中描述的主要流程(Fig. 4)。

虽然对本题而言，给定了 `N = 21`，但本组还是希望给出尽可能通用性高的实现。对于完整的 Shor 算法，其需要：

1. 简单判断给定的 `N` 是否是偶数，如是，显然 2 是一个非平凡的因子
2. 我们跳过下图中的步骤 2，在保证 `N` 为奇数的情况下，正式进行 Shor 算法的流程。首先需要随机选取一个候选的 `a`，其值范围在 `(1, N)` 开区间内。
3. 对于上述步骤选定的 `a`，显然有可能直接随机到目标的因子。因此需要简单判断 `a` 与 `N` 是否互质。如不是，则可以通过计算最大公因数得到一个非平凡的因子。例如针对 `N = 21, a = 6`，可以得到非平凡因子 `gcd(21, 6) = 3`。
4. 正式调用 period_finding 量子子程序求得目标 `a` 与 `N` 所对应的函数 `f(x) = a ^ x mod N` 的周期 `r`。这一部分才正式的与量子算法有关。
5. 根据求得的 `r` 计算 `gcd(a ^ (r / 2) + 1, N)` 与 `gcd(a ^ (r / 2) - 1, N)`，得到目标的因子。注意，按题目要求需排序后返回。


<br />
<div style="margin: 0 auto; width: 600px; text-align: center;">
    <img src="./images/4. shor.png" style="margin-bottom:15px"/>
    <label style="font-size: 17px;"><b>Fig. 4</b> Shor 算法主要流程</label>
</div><br />

则大体的框架如下：

```python
def solution() -> List[int]:
    N = 21
    res = shor_alg(N)

    return res

def shor_alg(N: int) -> List[int]:
    """ Shor algorithm using to factoring a composite number N

    """

    # Step 1. If N is even, 2 is the factor
    ...
    # Step 2. Randomly choose a number 1 < a < N
    a = random.randint(2, N - 1)
    ...
    # Step 3. Check the randomly choosed number a
    #
    # compute K = gcd(a, N), 
    # if K != 1, then K is thus a non-trivial factor of N
    # algorithm finished with returned [K, N / K]
    # 
    ...
    # Step 4. call quantum period-finding subroutine to find period r
    ...
    r = period_finding(a, N)
    # Step 5. calculate factors:
    factors = calculate_factors(a, N, r)

    return factors
```

其中的 `Shor_alg()` 的基本步骤按照上述文字描述的实现，基本是经典的处理，这里不做过多解释。 

下文会针对其中具体的 `period_finding()` 与 `calculate_factors()` 函数实现展开分析。


### 2.2 period_finding 的量子线路设计

### 2.2.1 period_finding 基本流程

period_finding 量子子程序需要构造类似 Fig. 1 所示的量子线路来求解函数 `f(x) = a ^ x mod N` 的周期 `r`。详细的步骤也有较多的参考资料，这里简单回顾下基本的流程。

<br />
<div style="margin: 0 auto; width: 500px; text-align: center;">
    <img src="./images/5. phase estimation.png" style="margin-bottom:15px"/>
    <label style="font-size: 17px;"><b>Fig. 5</b> 量子相位估计算法的一般线路结构</label>
</div><br />

首先，period_finding 子程序基本也可以看作是量子相位估计的一种应用，图 Fig. 5 是量子相位估计的一般线路结构，可以看出，大体结构上同 period_finding 是几乎相同的。其中的上半部的量子寄存器用于存储测量得到的相位估计的结果（这部分量子寄存器的个数就决定了相位估计结果的精度）。下半部分的辅助比特在 period_finding 语境下，用来存储 `f(x) = a ^ x mod N` 函数的运行结果。对于 `N = 21` 而言，显然不考虑特殊的优化的话，至少需要 5 个量子比特来表示。

用一个例子来解释 period_finding 求得函数 `f(x) = a ^ x mod N` 周期的基本过程。例如，对于 `a=7, N=15` 的情况，假设我们使用三个量子比特的精度来存储相位估计结果，并且测量得到了 110，即结果为二进制表示的 0.110，转为十进制为 0.75，可以写成精确的分数 3/4，则对应的分母就被认为是目标的周期 `r`。所以得到 `r = 4`。进而通过计算 `gcd(a ^ (r / 2) + 1, N)` 与 `gcd(a ^ (r / 2) - 1, N)` 可以得到非平凡的因子 5 和 3。

### 2.2.2 period_finding 线路设计

回到问题上来，如图 Fig. 1 所示，我们需要构造的量子线路大体上分为以下几个部分：

1. 准备一组量子寄存器用于存储相位估计的结果，准备另一组量子寄存器用于存储函数转换 `f(x) = a ^ x mod N` 的值。
2. 对第一组量子寄存器应用 Hadamard 变换制备叠加态。
3. 以第一组量子寄存器作为“输入”，以第二组量子寄存器作为函数转换 `f(x) = a ^ x mod N` 的输出，对所有量子比特运用一组变换（模指数线路部分）。
4. 对第一组量子寄存器进行 QFT 的逆变换，提取相位估计结果。

如前面所说的，其中的（3）步骤模指数线路的有效设计与实现是一个重要的瓶颈。相对较通用的实现可能会导致量子比特的使用数目偏多，门构造较复杂；而结合先验知识针对性设计的线路虽然可以极大的减少量子比特，但相应的也缺失了一定的可扩展性等。

此外，针对 `N = 21` 的情况，符合条件、需要构造对应量子线路的 `a`，亦即在 `(1, N)` 范围内与 `N` 互质的数包括：`[2,4,5,8,10,11,13,16,17,19,20]`，故如果要针对 `N = 21` 实现较完整的 Shor 算法，至少需要有能力设计与构建这几种 `a` 的取值情况下的模指数线路。

权衡后，本组采用了[<sup>3</sup>]所提供的一组针对 `N = 21` 的所有可能的 `a` 取值的模指数线路，见下图 Fig. 6。

<br />
<div style="margin: 0 auto; width: 1020px; text-align: center;">
    <img src="./images/6. transform_functions.png" style="margin-bottom:15px"/>
    <label style="font-size: 17px;"><b>Fig. 6</b> 对应 N=21 的各个可能的 a 取值的模指数线路</label>
</div><br />

这组线路的正确性此处不做解释，可以参考原文。需要指出的是，该组线路采用三个量子比特用作相位估计结果的第一组寄存器，使用五个量子比特存储函数变换的结果（如前面分析所说），共使用 8 个量子比特来完成分解 `N = 21` 的 Shor 算法的实现。

而受限于三个量子比特表征的相位估计结构精度有限，实验中对于诸如 `a = 5` 等几个个例运行效果不佳。后面我们采用重复实验等方式来减少由于这个问题导致的 `Shor_alg()` 调用不成功的情形发生。

### 2.2.3 程序实现

首先需要根据对应的 `N` 与 `a` 通过我们实现的 `create_program()` 函数返回构造好的对应线路。这里我们假设 `create_program()` 是一个有能力返回任意 `N` 与 `a` 组合的  period_finding 线路的函数，也没有写死 `N = 21` 以尽可能提升 API 的通用性。虽然事实上 `create_program()` 仅针对 `N = 21` 的情况返回对应 `a` 的线路（相应线路的模指数部分参考前述的[<sup>3</sup>]）。

组装线路部分主要就是按照运用 Hadamard 门制备叠加态 -> 应用对应的模指数转换函数线路 -> 逆 QFT 进行相位估计来构建量子线路。下图给出了 `N = 21, a = 19` 时构建的线路图为例：

<br />
<div style="margin: 0 auto; width: 1050px; text-align: center;">
    <img src="./images/10. shor_a=19.png" style="margin-bottom:15px"/>
    <label style="font-size: 17px;"><b>Fig. 7</b> 对应 N=21, a=19 构建的线路图</label>
</div><br />

得到线路后，将运行并测量以得到相位估计结果，即完成前述例子中的 110 => 0.75 的部分。此部分尤其需要注意读取结果的顺序（下图中 Step 2）。

```python
def period_finding(a: int, N: int) -> int:
    """ period_finding subroutine called in shor-algorithm
    
    """
    # Step 1. init period-finding QProg corresponding to a & N
    #
    # while actually in this experiment just implement the case N=21
    # 
    prog = create_program(a, N)

    # Step 2. Measure it (first 3 phase_bits_qubits qubits) => get phase
    #
    # attention the reading order: qubits[i] => cbits[phase_bits_qubits - i - 1]
    # 
    for i in range(phase_bits_qubits):
        prog << Measure(qubits[i], cbits[phase_bits_qubits - i - 1])

    result = run_with_configuration(prog, cbits, 1)
    print(f'  result: {result}') # like {"101": 1}

    # Convert the reading bits to phase
    #
    # eg. {"101": 1} => ["101"] => 5 => 5 / (2^3) = 0.625
    # 
    phase = int(list(result)[0], 2) / (2 ** phase_bits_qubits)
    print(f'   - corresponding phase: {phase}')

    # Step 3. calculate period r (by continued fraction expansion)
```

完成上面步骤后，将对得到的小数值（相位估计结果）进行进一步的处理以提取周期 `r`。这步重要的是需要将求得的小数转换成 `s/r` 的分数形式，通常采用连分式计算的方式来逼近的得到近似的分数表示，进而得到周期 `r`。Fig. 8 给出了一个连分式逼近的比较直观的例子。

<br />
<div style="margin: 0 auto; width: 900px; text-align: center;">
    <img src="./images/7. continued fraction.png" style="margin-bottom:15px"/>
    <label style="font-size: 17px;"><b>Fig. 8</b> 连分式算法举例</label>
</div><br />

我们可以借助 python 已有的标准库模块 `fractions` 中的 `Fraction.limit_denominator()` 函数来实现。事实上背后的原理也是利用连分式算法。

```python
    ...
    # Step 3. calculate period r (by continued fraction expansion)
    # 
    frac = Fraction(phase).limit_denominator(N)
    r = frac.denominator
```

### 2.2.4 周期 `r` 的校验

到上面的步骤似乎已经能够得到周期 `r` 了，但注意，根据 Shor 算法步骤的要求（见图 Fig. 9[<sup>5</sup>]），对求得的周期 `r` 有相关的限制，见图中框出的部分。

<br />
<div style="margin: 0 auto; width: 900px; text-align: center;">
    <img src="./images/8. validate.png" style="margin-bottom:15px"/>
    <label style="font-size: 17px;"><b>Fig. 9</b> shor 算法步骤 - 对求得的周期 r 的验证</label>
</div><br />

针对 Shor 算法分解 `N = 21` 的实验与讨论[<sup>1</sup>] [<sup>6</sup>]都指出，对于 `a = 4` 等情况下，调用能返回正确解的场合，对应的 `f(x) = a ^ x mod N` 周期为 `r = 3`，见 Fig. 10。然而如图 Fig. 9 中蓝色框标明的一样，Shor 算法的原本设计描述中指出期待求得的周期 `r` 为偶数，事实上，要求如果求得 `r` 为奇数，意味着调用失败，需要回退到第一步。

<div style="margin: 0 auto; width: 470px;
text-align: center;">
    <img src="./images/9. periods.png" style="margin-bottom:15px; position: relative; left: 40px"/>
    <label style="font-size: 17px;"><b>Fig. 10</b> N=21 时 f(x) 的周期情况</label>
</div><br /><br />

针对上述问题，在 `r = 3` 结论正确的前提下，显然有必要放宽 `period_finding` 过程中针对 `r` 的限制。事实上，[<sup>1</sup>] [<sup>6</sup>] [<sup>7</sup>] 都指出了，在一定条件下，`r` 为奇数并不影响经典 `period_finding` 过程的求解。[<sup>6</sup>] 给出了一个简单的筛选方法，即当选取的 `a` 本身为完全平方数时，即使 `r` 为奇数也并不影响后续的因子计算。主要原因在于当 `a` 为完全平方数时，可以保证即使 `r` 为奇数，`a ^ (r / 2)` 也是整数。从而不妨碍 Fig. 9 所示的第 7 步的正常计算。

在这个条件下，Fig. 10 中所示 `a = 4` 与 `a = 16` 的情况(均有 `r = 3`)可以正常的计算。

为此，在 `period_finding` 的实现中对合法求得的 `r` 为奇数的情况不做排斥与回退。但依旧需要检查 `a ^ r mod N == 1` 的条件。另外，由于在 `a = 4` 为例的情况下，连分式得到的结果(eg. `r = 8`) 可能并不满足，所以简单的处理方式可以保留连分式过程中的每个子结果，记录每个可能的 `r` 并逐个判断[<sup>6</sup>]。例如对于求得 `phase=5/8, a=4`，连分式的逼近过程可以得到 `{0, 1, 1/2, 2/3, 5/8}` 则可以分别判断 `r=8, r=3, ...`。

从实现上，借助 `Fraction().limit_denominator()` 来求得连分式结果，通过不断缩小 `limit_denominator` 的参数，得到逼近过程中的每个解。并逐次判断是否合法。

按上述讨论的，逼近式的求解合法的 `r`，这一步主要验证了 Fig. 9 中红框所示部分。

```python
    ...
    # Step 3. calculate period r (by continued fraction expansion)
    # 
    # eg. a = 4, phi = 0.625 gives the convergents: {0, 1, 1/2, 2/3, 5/8}
    # while r = 8 is invalid (4^8 mod 21 = 16 !== 1)
    # actually in this case only r = 3 is valid (4^3 mod 21 == 1)
    # 
    # according to the above discussion we'll narrow the limit unitl getting a valid r
    # or narrow to limit = 0 which means fail to find period. 
    #
    limit = N

    while True:
        frac = Fraction(phase).limit_denominator(limit)
        r = frac.denominator

        # simply check period
        if (a ** r) % N == 1:
            break
        else:
            # narrow limit to calculate new period
            limit = r - 1
            if limit <= 0:
                print(f'\n  Rewrite to fraction: {frac}, find period failed')

                return
```

对于 Fig.9 中蓝框部分的验证，如前面所讨论的，当 `r` 为奇数，但 `a` 本身为完全平方数时，被认为是合法的，不做排斥。

```python
    ...
    # re-check calculated r
    # 
    # a itself is a perfect square thus Shor still works
    # cite: https://arxiv.org/pdf/2103.13855.pdf
    # 
    if (r % 2 != 0 and int(math.sqrt(a)) ** 2 != a) or \
        (int(a ** (r / 2)) % N == -1):
        print(f'\n  Rewrite to fraction: {frac}, find period failed')

        return
```

### 2.3 calculate_factors() 函数实现

到这步，主要就是基于 Fig. 9 所示的第 7 步根据已求得的周期 `r`，计算因子。当然，计算得到的结果可能并不是非平凡的解，因此需要再经过一次检验，过滤掉平凡的 `[1, 21]`。最后，根据接口要求，需要返回排序后的结果。

```python
def calculate_factors(a: int, N: int, r: int) -> List[int]:
    """ calculate factors based on calculated period r.

    """
    # According to Shor algorithm, calculate gcd(a^(r/2) - 1, N) and gcd(a^(r/2) + 1, N)
    # 
    guesses = [math.gcd(int(a ** (r / 2)) - 1, N), math.gcd(int(a ** (r / 2)) + 1, N)]
    
    print(f'  calculate final guesses: {guesses}')
    
    # need to check the calculated guesses numbers.
    # 
    ...
    
    return sorted(factors)

```

## 3. 优化

### 3.1 重复实验提升成功率

前面也提到，对于部分 `a` 的取值，在当前的三个量子比特来表征量子相位估计结果的精度的条件下可能会出现失败，加上量子测量本身存在概率性，取单次测量结果来进行计算也有一定的概率失败，故从提升程序成功率的角度来看，重复实验是最简单的改进办法。

包括在失败的情况下重复调用 `shor_alg`：

```python
    ...
    # Cause the shor algorithm might fail => attempt some rounds
    MAX_ROUND = 8
    
    N = 21

    for round in range(MAX_ROUND):
        print(f'Attempt call shor-algorithm: round {round + 1}\n')

        res = shor_alg(N, ctx = ctx)
        ...
```

以及在 period_finding 子程序中可以重复多轮求取周期，以提升成功率。为此，适当改造了 `period_finding()` 和 `calculate_factors()` 函数，返回一个变量来判断调用成功与否。

```python
    # Step 3. call quantum period-finding subroutine to find period r
    #
    # because each time the running result will affected by nondeterministic measurements
    # (during period-finding, will measure the phase result to calculate period) 
    # for selected {a}, will try {MAX_ATTEMPT} to calculate period repeatly
    # 
    MAX_ATTEMPT = 20
    attempt = 1
    
    while True:
        print(f'Execute shor algorithm - attempt time: [{attempt}]')

        # call period-finding subroutine
        valid, r = period_finding(a, N, ctx)
        if valid:
            # valid period, and then calculate factors:
            ok, factors = calculate_factors(a, N, r)
            if ok:
                ...
            
        attempt += 1
        ...
```

在上述简单的修改之后，整体的成功率有了明显的提升。

### 3.2 Toffoli 门近似

作为优化的一环，完整的 Toffoli 门的实现需要 15 个量子逻辑门，可以使用带有相对相移的近似 Toffoli 门[<sup>8</sup>]实现以减少逻辑门的使用个数(见下图 Fig. 12)。已有的实验[<sup>3</sup>]表明此种近似的 Toffoli 门可应用于 Shor 算法的实现。

<div style="margin: 0 auto; width: 770px;
text-align: center;">
    <img src="./images/12. approximate.png" style="margin-bottom:15px; position: relative; left: 40px"/>
    <label style="font-size: 17px;"><b>Fig. 11</b> Toffoli 门近似</label>
</div><br /><br />

在本地进行测试，针对不展开 Toffoli 门、展开 Toffoli 门、采用相对相移近似 Toffoli 门分别运行 5000 次，并统计尝试次数(每次随机选择一个 `a` 记一次)，平均尝试次数与成功率。

<div style="margin: 0 auto; width: 770px;
text-align: center;">
    <img src="./images/13. table.png" style="margin-bottom:15px; position: relative; left: 40px"/>
    <label style="font-size: 17px;"><b>Fig. 12</b> 测试结果</label>
</div><br /><br />

实验表明，该 Shor 算法的实现对于分解特定的问题 `N = 21` 有较好的表现。此外，实现结果也没有明显表现出选择或不选择展开 Toffoli 门对实验(尝试次数、平均尝试次数、成功率)的影响。

可以认为采用相对相移版本的 Toffoli 门近似可以在近乎相同的表现下，节省量子门数量。故最终解选择采用上述相对相移的 Toffoli 门近似方案。

## **源码说明**

所给的源码在 python == 3.7.0 环境下可以运行。下面给出了本组实现的结果（未使用 Toffoli 近似时的结果）：

<div style="margin: 0 auto; width: 1100px;
text-align: center;">
    <img src="./images/14. all.png" style="margin-bottom:15px; position: relative; left: 40px"/>
    <label style="font-size: 17px;"><b>Fig. 13</b> 针对 N=21 的各个 a 取值的模指数线路</label>
</div><br /><br />





[<sup>1</sup>]:https://arxiv.org/pdf/1310.6446v2.pdf
[<sup>2</sup>]:http://mmrc.amss.cas.cn/tlb/201702/W020170224608149940643.pdf
[<sup>3</sup>]:https://arxiv.org/pdf/1202.6614.pdf
[<sup>4</sup>]:https://www.scottaaronson.com/qclec/21.pdf
[<sup>5</sup>]:https://en.wikipedia.org/wiki/Shor%27s_algorithm
[<sup>6</sup>]:https://arxiv.org/pdf/2103.13855.pdf
[<sup>7</sup>]:https://www.nature.com/articles/nphoton.2012.259
[<sup>8</sup>]:https://arxiv.org/pdf/1508.03273.pdf
