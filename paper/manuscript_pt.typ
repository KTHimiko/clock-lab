// Versão em português. Regra tipográfica seguida aqui: o modo matemático do
// Typst trata a vírgula como separador e insere espaço depois dela, de modo que
// $0,907$ sairia como "0, 907". Com decimal em vírgula isso apareceria no
// documento inteiro. Portanto: símbolos e fórmulas em matemática, números
// soltos em texto corrido.
#set page(paper: "a4", margin: (x: 2.3cm, y: 2.4cm), numbering: "1")
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "pt")
#set par(justify: true, leading: 0.62em, first-line-indent: 0pt, spacing: 0.9em)
#set heading(numbering: "1.1")
#show heading: it => block(above: 1.4em, below: 0.7em)[
  #set text(size: if it.level == 1 { 12pt } else { 10.5pt }, weight: "bold")
  #if it.numbering != none [#counter(heading).display(it.numbering)#h(0.6em)]
  #it.body
]
#show figure.caption: it => [
  #set text(size: 9pt)
  #set par(justify: true)
  *#it.supplement #context it.counter.display(it.numbering).* #it.body
]
#set figure(gap: 0.9em, supplement: [Figura])

#align(center)[
  #block(text(size: 15pt, weight: "bold")[
    A correção de composição celular na idade epigenética\
    pode injetar o confundimento que deveria remover
  ])
  #v(0.2em)
  #block(text(size: 11.5pt)[
    Um índice de risco de transporte, e um conserto de uma linha
  ])
  #v(1.1em)
  #text(size: 10.5pt)[Luan Ivepe]
  #v(0.2em)
  #text(size: 9.5pt, style: "italic")[Pesquisador independente]
  #v(0.2em)
  #text(size: 9pt)[#link("mailto:luanivepe@gmail.com")[luanivepe\@gmail.com]]
  #v(0.2em)
  #text(size: 9pt)[Rascunho de preprint — #datetime.today().display("[day]/[month]/[year]")]
  #v(0.5em)
  #text(size: 8.5pt, style: "italic")[
    Versão em português. A versão de referência, para submissão, é a inglesa.
  ]
]

#v(1.0em)

#block(inset: (x: 1.2em), [
  #text(weight: "bold")[Resumo] #h(0.6em)
  A aceleração de idade epigenética em sangue é rotineiramente corrigida pela
  composição de células imunes, e essa correção é quase sempre estimada e
  aplicada dentro da mesma coorte. Mostramos que, quando os coeficientes dessa
  correção são estimados numa coorte e aplicados em outra — a situação sempre
  que uma correção publicada é reutilizada — o procedimento pode *adicionar*
  sinal de composição em vez de removê-lo. Ajustando em quarenta amostras e
  transportando, o termo de composição residual sobe 12,2 pontos da variância de
  aceleração de idade e a correção fica pior do que não corrigir em 89% dos
  sorteios; já coeficientes de composição sem informação alguma ficam em +0,5%,
  com intervalo interquartil atravessando o zero. O dano vem de coeficientes
  reais mal estimados, não de ruído. A falha é invisível dentro da coorte de
  ajuste, onde a correção zera o termo de composição por construção. Nem o
  tamanho amostral nem o condicionamento da matriz de composição da coorte de
  ajuste explicam o fenômeno — a coorte que falha pior é a mais bem condicionada
  do conjunto. O que explica é o termo padrão de excesso de risco de mínimos
  quadrados sob _covariate shift_, $(sigma^2 \/ n) dot tr(Sigma_"ajuste"^(-1)
  Sigma_"teste")$, que chamamos de índice de transporte. Ele ordena 63
  configurações direcionadas em quatro coortes de ajuste com $rho$ de Spearman
  de 0,907, é assimétrico onde distâncias simétricas não são, e prevê a pior
  falha observada com meio ponto de erro. Penalizar os coeficientes de
  composição por ridge limita a amplificação diretamente e mantém os doze pares
  direcionados em zero ou abaixo, inclusive um transporte que custa +15,9 pontos
  sem penalidade. A validação cruzada comum na coorte de ajuste *não* dimensiona
  a penalidade corretamente justamente na coorte que falha. Reportamos o índice
  como diagnóstico computável de antemão e a penalidade como padrão.
])

= Introdução

Sangue é uma mistura, e a mistura muda com a idade. Linfócitos caem, células
mieloides sobem, e compartimentos naive dão lugar aos de memória. Como cada
subtipo de leucócito carrega o próprio metiloma, um relógio epigenético aplicado
a sangue total lê em parte a idade das células e em parte o censo do tubo
@jaffe2014.

O tamanho dessa contaminação já está bem caracterizado. Em resolução de doze
tipos, a composição imune explica de 13% a 34% da variância de aceleração de
idade, dependendo do relógio @zhang2024, e células T CD8 naive leem 15 a 20 anos
mais novas que CD8 de memória efetora do mesmo doador @tomusiak2024. Existem
duas respostas. Uma projeta o confundimento para fora já na seleção de CpGs
@tomusiak2024. A outra — muito mais comum, e o objeto deste artigo — subtrai o
efeito depois: regride a idade do relógio sobre a idade cronológica e as
proporções celulares estimadas, e guarda o resíduo. Isso é a aceleração de idade
epigenética intrínseca, e é computada dentro de qualquer conjunto de dados que
esteja à mão.

É esse último detalhe que importa. *Uma correção estimada e avaliada na mesma
coorte não pode ser validada ali.* O resíduo de um ajuste de mínimos quadrados é
ortogonal aos seus preditores por construção, então o termo de composição é
levado a exatamente zero independentemente de quanto os coeficientes valham.
Qualquer avaliação de se a correção funciona precisa portanto ser feita fora da
coorte — e o próprio artigo de referência da área, que compara oito métodos de
correção de tipo celular, não testa transferência entre conjuntos de dados para
nenhum deles #footnote[Genome Biology 17:84 (2016), "An evaluation of methods
correcting for cell-type heterogeneity in DNA methylation studies"; PMC4855979.
Lista de autores a conferir antes da submissão.].

Transferência não é hipótese. Sempre que uma correção publicada é reutilizada,
que uma coorte clínica de algumas dezenas é corrigida com coeficientes de um
estudo maior, ou que dois grupos comparam estimativas corrigidas, coeficientes
atravessam uma fronteira de coorte. Perguntamos o que acontece quando
atravessam.

= Métodos

== Coortes, painéis de referência e relógios

Quatro séries públicas de sangue total em 450k com idade cronológica: GSE40279
($n = 656$, idades 19–101), GSE61151 ($n = 184$), GSE50660 ($n = 464$) e
GSE42861 ($n = 689$). As proporções celulares foram estimadas por mínimos
quadrados não negativos com restrição, contra dois painéis de referência
construídos aqui: um painel de seis tipos a partir de Reinius et al. (GSE35069)
e um de doze tipos a partir da referência FlowSorted.BloodExtended.EPIC
(GSE167998) @salas2022. Nosso painel de doze tipos usa uma seleção de sondas por
estatística t mais simples que a biblioteca IDOL publicada, e recupera as
proporções conhecidas das doze misturas reconstruídas com $r$ de 0,79 e erro
absoluto médio de 0,027; todo número de composição aqui deve ser lido como um
piso do que a biblioteca publicada alcança. O canal de monócito do painel corre
alto, e nenhum coeficiente de monócito por tipo é interpretado.

Três relógios com coeficientes publicados e abertos sustentam a análise: Horvath
2013, Levine 2018 (PhenoAge) e Horvath 2018. Hannum 2013 foi treinado no GSE40279
e é excluído sempre que essa coorte aparece.

== A correção, e como ela é pontuada

A correção é a padrão. Na coorte de ajuste, ajustamos
$y = beta_0 + beta_1 a + C beta_c$, em que $y$ é a idade do relógio, $a$ a idade
cronológica e $C$ as proporções estimadas (uma coluna descartada, já que somam
um); guardamos $beta_c$ e aplicamos
$y - (C_"teste" - macron(C)_"ajuste") beta_c$ na coorte de teste. A idade
cronológica permanece no modelo de ajuste — composição deriva com a idade, e
omiti-la arrancaria envelhecimento real junto com o hemograma — mas só os
coeficientes de composição viajam, então aplicar a correção não exige idade.

O desfecho é o termo de composição que sobra no resíduo de idade da coorte de
teste: o incremento em $R^2$ ao acrescentar composição a um modelo da idade do
relógio sobre a idade cronológica, menos um nulo de permutação, como fração da
variância de aceleração de idade não corrigida. *Ajuste e medição usam painéis
diferentes* — ajustado com doze tipos, medido com seis — porque pontuar uma
correção com a própria representação dela é pontuá-la contra o próprio relato do
que removeu. Reportamos $Delta$ = depois $-$ antes; positivo significa que a
correção deixou o relógio pior do que não corrigir.

== O índice de transporte

Escreva o erro de estimação como $e = hat(beta)_c - beta_c$. A correção
transportada deixa $-C_"teste" e$, cuja variância na coorte de teste é
$e' Sigma_"teste" e$. Para mínimos quadrados,
$"Cov"(e) = sigma^2 \/ n dot Sigma_"ajuste"^(-1)$, logo

$ EE["dano"] prop (sigma^2 / n) dot tr(Sigma_"ajuste"^(-1) Sigma_"teste") . $

Chamamos $tr(Sigma_"ajuste"^(-1) Sigma_"teste") \/ n$ de *índice de transporte*.
Ambas as matrizes de segundo momento são tomadas depois de retirar intercepto e
idade cronológica, que é o espaço onde os coeficientes vivem.

*Não reivindicamos novidade para essa quantidade.* Ela é o termo padrão de
excesso de risco de mínimos quadrados sob _covariate shift_ — a situação em que
a distribuição das covariáveis muda entre treino e aplicação — e o modo de falha
que ela descreve tem nome nessa literatura: _spectral inflation_, direções que
carregam pouca variação no treino e mais na avaliação @spectral2023. Há também
teoria sobre escolher a regularização ridge sob _covariate shift_ @patil2024. O
que é novo aqui é a observação de que uma correção epidemiológica de uso corrente
está exposta a isso, e a medição de quando o problema morde.

Para ridge com penalidade $lambda$, o termo correspondente é
$tr(M Sigma_"ajuste" M Sigma_"teste") \/ n$, com
$M = (Sigma_"ajuste" + lambda I)^(-1)$. As penalidades são citadas como $alpha$,
um multiplicador livre de escala do autovalor médio do produto cruzado
padronizado, de modo que o mesmo $alpha$ signifique a mesma coisa em todo $n$ e
em toda coorte. Os coeficientes são ajustados por Frisch–Waugh–Lovell sobre
colunas padronizadas e devolvidos à escala original; em $alpha = 0$ isso
reproduz exatamente o ajuste sem penalidade.

== Subamostragem

Onde uma coorte de ajuste é subamostrada, os sorteios são *estratificados por
decil de idade*. Um sorteio não estratificado de quarenta, numa coorte que vai de
19 a 101 anos, estreita e desloca a faixa etária, o que confundiria tamanho
amostral com faixa de treino. Ao longo da grade, a idade mediana de um sorteio
fica a menos de 0,5 ano da coorte cheia, e a amplitude entre os percentis 10 e 90
dentro de 2% dela.

= Resultados

== Uma correção transportada pode ser pior do que não corrigir

#figure(
  image("figures/pt/fig1_curva_n.png", width: 95%),
  caption: [*A correção é nociva abaixo de aproximadamente 160 amostras de
  ajuste.* Composição que sobra no resíduo de idade depois do transporte, como
  fração da variância de aceleração de idade; acima de zero, a correção deixou o
  relógio pior do que não mexer nele. Correções ajustadas em subamostras
  estratificadas por idade do GSE40279 e levadas a três coortes externas, 30
  sorteios por ponto, três relógios. As linhas cinza são as três coortes de teste
  separadamente. A régua laranja é o mesmo transporte com coeficientes de
  composição ajustados sobre composição embaralhada, que não carregam informação
  por construção.],
) <fig1>

Ajustando a correção em subamostras do GSE40279 e transportando-a para três
coortes externas (@fig1), a mediana de $Delta$ cruza zero entre $n$ de 160 e
$n$ de 184. Abaixo disso a correção padrão é pior do que não fazer nada, e no
fundo da faixa não é por pouco: com quarenta amostras, $Delta$ chega a +12,2% e a
correção é nociva em 89% dos sorteios. Expressa contra o sinal de composição que
estava presente, a correção em $n$ igual a 40 não remove 0% dele — ela *adiciona
2,7 vezes* o que havia. Com as 656 completas, remove 61% fora da coorte.

A linha de referência importa para a interpretação. Coeficientes de composição
ajustados sobre composição *embaralhada* não carregam informação, e
transportá-los produz uma mediana de $Delta$ de +0,5%, com intervalo
interquartil de −0,8% a +2,1% sobre 180 sorteios — uma distribuição que atravessa
o zero. Coeficientes sem sentido são praticamente inofensivos. São os
coeficientes reais, estimados em amostras de menos, que causam o dano, e causam
cerca de vinte vezes mais dele do que o ruído causa.

A mediana não é a história toda. Tomando o espalhamento entre sorteios *dentro*
de cada célula relógio-por-coorte, o quartil superior da célula mediana só cai
abaixo de zero por volta de $n approx 320$: o sorteio típico começa a ajudar
perto de 160 a 184, mas um sorteio azarado continua prejudicando bem depois
disso.

== Nem tamanho amostral nem condicionamento explicam

Tamanho amostral é causa forte e não suficiente. Ajustar em 184 amostras do
GSE40279 e transportar dá $Delta$ de −0,9% — levemente benéfico — enquanto
ajustar nas 184 do próprio GSE61151 e transportar para o GSE40279 dá +5,0%,
com os três relógios prejudicados. Alguma coisa na coorte de ajuste, além do
tamanho dela, está fazendo o trabalho.

O candidato natural é o condicionamento: proporções celulares colineares dão
coeficientes instáveis, e fatores de inflação de variância acima de 100 são
documentados como capazes de inverter o sinal da maioria dos coeficientes nesse
contexto @meredith2019. No nosso painel de doze tipos, o maior fator de inflação
de variância vai de 215 a 302, pior que os 113,7 relatados com seis tipos. Mas o
condicionamento mal se move com $n$ — o número de condição do bloco de composição
residualizado é 65,9 em $n$ igual a 40 e 54,5 em $n$ igual a 656 — e, dentro de
cada célula $n$-por-relógio-por-coorte, não carrega informação alguma sobre o
dano ($rho$ de Spearman mediano de −0,025, positivo em 48% das células, um cara
ou coroa).

Decisivamente: *o GSE61151, a coorte cuja correção transportada falha pior, é a
mais bem condicionada do conjunto* (número de condição 40,8, maior VIF 117,8,
contra 54,4 a 65,9 e 212 a 302 para toda subamostra do GSE40279, inclusive a
coorte cheia). O condicionamento da coorte de ajuste sozinho não pode ser o
mecanismo.

== O índice de transporte ordena o dano

#figure(
  image("figures/pt/fig2_indice.png", width: 95%),
  caption: [*O índice de transporte ordena 63 configurações direcionadas.* Cada
  ponto é uma configuração (coorte de ajuste, coorte de teste, tamanho de
  ajuste); a forma do marcador dá a coorte de ajuste. Spearman de 0,907.
  Nenhuma configuração abaixo de um índice de 0,05 foi nociva; 81% acima dele
  foram.],
) <fig2>

O índice é uma propriedade conjunta das duas coortes, que é exatamente o motivo
de medir $Sigma_"ajuste"$ sozinha não ter dado nada: uma $Sigma_"ajuste"$ bem
condicionada pode ainda assim ter $Sigma_"ajuste"^(-1)$ amplificando justamente
as direções em que $Sigma_"teste"$ guarda variância.

Em 63 configurações direcionadas montadas a partir de quatro coortes de ajuste, o
índice ordena o dano mediano com $rho$ de 0,907 ($p = 2 times 10^(-24)$; @fig2),
e faz isso dentro de cada coorte de ajuste separadamente — 0,858, 0,951, 0,907 e
0,897 — de modo que a ordenação não é propriedade da escada de subamostras de uma
coorte. O dano cruzou zero perto de um índice de 0,05: *nenhuma* das 21
configurações abaixo desse valor foi nociva, e 81% das acima foram. Entre 0,05 e
0,16 ocorrem os dois desfechos, então o limiar é um piso de um lado só, não um
ponto de travessia.

Um resultado negativo pertence a esta seção. O índice *não* prevê o dano de um
sorteio individual: com $n$ e coorte de teste fixos, o $rho$ mediano é 0,021 em
56% das células, pior que a grandeza que ele deveria superar. A teoria prevê
isso. Com $n$ e a coorte de teste fixos, $Sigma_"ajuste"$ é quase idêntica entre
sorteios, então o índice mal varia, enquanto o $e' Sigma_"teste" e$ realizado é
do tipo qui-quadrado com onze graus de liberdade. Um preditor quase constante não
consegue acompanhar um desfecho dominado por ruído de realização. O índice prevê
dano esperado entre configurações, não o sorteio à sua frente.

#figure(
  image("figures/pt/fig3_assimetria.png", width: 95%),
  caption: [*O índice é assimétrico, e o dano também.* Cada par reversível
  ajustado nas duas direções em $n = min(n_A, n_B)$, de modo que só a covariância
  de ajuste muda. O laranja marca a direção que o índice chama de pior; ela está
  à direita do marcador vazado em cinco de seis. No sexto
  (GSE40279 / GSE42861) os dois índices são 0,0275 e 0,0271 — o índice previu
  indiferença e foi pontuado como se tivesse previsto uma direção.],
) <fig3>

Distâncias simétricas entre estruturas de covariância de composição não dão conta
desse trabalho, e a assimetria é a razão. Ajustado em $n$ casado, de modo que só
a covariância de ajuste difere, o índice nomeia a direção pior em cinco de seis
pares reversíveis (@fig3; binomial $p$ de 0,109, que não é significância e não é
apresentado como tal — com seis pares, nenhum limiar simultaneamente passa de
0,05 e tolera um erro). O erro é um par cujos dois índices são 0,0275 e 0,0271: o
índice previu indiferença e foi pontuado como se tivesse previsto direção.

Para a configuração que derrotou nossas tentativas anteriores, o índice de
GSE61151 $arrow$ GSE40279 é 0,174 contra 0,038 do seu próprio reverso, uma razão
de 4,6. Um índice perto de 0,18 corresponde a cerca de +4,6% de dano na
@fig2; o valor medido é +5,0%.

Uma simulação fecha o mecanismo. Sorteando duas coortes sintéticas com o *mesmo*
vetor de coeficientes verdadeiros, o mesmo ruído e o mesmo $n$, de modo que a
única coisa que pode dar errado é estimação atravessando um desencontro de
covariância, trocar qual covariância faz o ajuste move o resíduo de +21,2% para
+84,6% em $n$ igual a 40. Os resíduos observados regridem sobre os previstos com
inclinação 0,92 pela origem e $r$ de 0,815.

== Penalizar os coeficientes remove a falha

#figure(
  image("figures/pt/fig4_conserto.png", width: 95%),
  caption: [*Ridge mantém os doze transportes direcionados em zero ou abaixo.*
  Cada linha é um par direcionado de coortes em $n = min(n_A, n_B)$; a seta vai
  do ajuste sem penalidade até $alpha = 3$.],
) <fig4>

Ridge substitui $Sigma_"ajuste"^(-1)$ por $(Sigma_"ajuste" + lambda I)^(-1)$ e
limita a amplificação diretamente. Uma penalidade só conta como conserto se
remove o dano *e* preserva o benefício: encolher uma correção até zero remove o
dano e não é método. Exigindo, antes da rodada, que a penalidade não seja pior
que não corrigir em $n$ igual a 40 *e* retenha ao menos 70% do benefício sem
penalidade em $n$ igual a 656, exatamente um valor da nossa grade qualifica,
$alpha = 3$. A curva dele é plana — entre −1,4% e −2,0% ao longo de uma faixa de
dezesseis vezes em tamanho de ajuste — e é esse o ponto: a correção deixa de
depender de quão grande era a coorte de ajuste. Já $alpha = 10$ é a esquina
degenerada que o critério foi escrito para pegar: norma dos coeficientes 13, dano
preso perto de −0,9% em toda parte, inofensivo e inútil.

O mesmo $alpha$, obtido independentemente a partir da coorte real que falhou, leva
aquele transporte de +5,0% para −0,5%. Nos doze pares direcionados em $n$
casado, mantém a mediana em zero ou abaixo (@fig4), inclusive
GSE61151 $arrow$ GSE42861, o pior transporte que encontramos, que vai de +15,9%
para −2,5%. Não é de graça: onde a correção sem penalidade já transporta bem, a
penalidade custa parte do benefício — GSE40279 $arrow$ GSE42861 abre mão da maior
parte dele, de −2,8% para −0,7%. Penalizar é seguro, e seguro tem prêmio nos
casos que não precisavam.

*A validação cruzada não basta para dimensioná-la.* Leave-one-out na coorte de
ajuste escolhe sensatamente em subamostras do GSE40279 e, na coorte que de fato
quebra o transporte, escolhe $alpha = 0.3$ e deixa +2,9% do dano de pé. A
validação cruzada otimiza predição dentro da coorte de ajuste, e esse objetivo
não sabe que os coeficientes estão prestes a ser embarcados para outro lugar.

== O diagnóstico sobrevive a duas equipes e dois painéis de referência

Tudo acima usa uma mesma deconvolução nas duas pontas de cada transporte, o que
não é a situação de que o resultado trata. Estimando a composição da coorte de
ajuste a partir de uma referência e a da coorte de teste a partir de outra — os
mesmos seis rótulos, doadores, sondas e geração de array diferentes — os dois
painéis concordam bem ($r$ mediano por tipo de 0,936), e o índice ainda ordena
as configurações com $rho$ de 0,776, contra 0,793 para painéis casados.
Duas referências discordando sobre um tipo celular *são* um desencontro de
covariância, então o índice enxerga sem precisar ser avisado. Ridge com
$alpha = 3$ mantém os doze pares em zero ou abaixo também no braço desencontrado.

= Discussão

Três coisas decorrem disso para a prática.

*Primeiro, uma correção estimada dentro da coorte não pode ser validada nessa
coorte.* Isso é aritmética, não achado empírico, e se aplica a todo artigo que
reporta uma aceleração de idade intrínseca sem uma checagem fora da coorte. A
perfeição aparente da correção em casa é garantida e não carrega informação.

*Segundo, coortes de ajuste pequenas não são meramente ruidosas.* Um vetor de
coeficientes sem informação nenhuma praticamente não faz nada quando
transportado. Coeficientes reais estimados em quarenta amostras adicionam 2,7
vezes o sinal de composição que estava presente. A falha não é uma versão
diluída do sucesso.

*Terceiro, o risco é computável de antemão.* O índice de transporte precisa
apenas das proporções estimadas das duas coortes e do tamanho de ajuste — nenhum
desfecho, nenhuma segunda coorte de validação, nenhum rótulo retido. Abaixo de
0,05 não observamos dano em 63 configurações. Não leríamos isso como garantia, e
sim como um limiar que vale calcular antes de reutilizar uma correção publicada.

Nossa recomendação é mais estreita que "use ridge e valide cruzadamente".
Penalize — e não deixe a validação cruzada na coorte de ajuste decidir sozinha
quanto.

= Limitações

Quatro coortes, todas de sangue total em 450k, todas adultas, e doze pares
direcionados. A contagem de assimetria é seis. Os dois painéis de referência são
construídos aqui com seleção de sondas mais crua que as bibliotecas publicadas,
então a discordância entre eles é um piso do que dois laboratórios reais
produziriam. O valor $alpha = 3$ é este painel com estes relógios; o que
esperaríamos que transportasse é "penalize", não o número. O índice é derivado
para mínimos quadrados com um termo de composição linear corretamente
especificado, e a simulação que confirma sua constante assume exatamente isso. E
o índice prevê dano esperado entre configurações, não o sorteio individual — um
resultado negativo que reportamos em vez de enterrar.

Nada aqui diz respeito a se relógios epigenéticos medem algo real. Esta é uma
afirmação sobre uma correção aplicada a eles.

= Disponibilidade de dados e código

As cinco séries são públicas: GSE40279, GSE61151, GSE50660, GSE42861 (sangue
total), GSE35069 e GSE167998 (painéis de referência). Nenhum dado de terceiros é
redistribuído. O código de análise, o registro etapa a etapa completo — incluindo
toda conclusão posteriormente derrubada — e os scripts das figuras estão em
#link("https://github.com/KTHimiko/clock-lab")[github.com/KTHimiko/clock-lab].

#v(0.8em)
#block(fill: luma(245), inset: 9pt, radius: 3pt, width: 100%, [
  #text(size: 9pt, weight: "bold")[Antes da submissão — a conferir]
  #text(size: 9pt)[
    #list(
      spacing: 0.45em,
      [A lista de autores do artigo de benchmark de oito métodos da Genome
       Biology 2016 (PMC4855979) e a do preprint de _covariate shift_
       arXiv:2312.17463 não foram confirmadas, e estão marcadas como tal em vez
       de adivinhadas.],
      [O repositório está privado no momento; precisa estar público antes de o
       preprint citá-lo.],
      [Falta uma declaração de reprodutibilidade e as versões exatas dos
       pacotes.],
    )
  ]
])

#bibliography("refs.bib", title: "Referências", style: "nature")
