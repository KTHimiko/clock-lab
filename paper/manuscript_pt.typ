// Versao em portugues. Regra tipografica: no modo matematico do Typst a virgula
// ganha espaco depois ("0, 907"), entao numeros decimais ficam em texto corrido
// e so simbolos e formulas ficam em matematica.
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
    Correções de composição celular da idade epigenética podem,\
    quando transportadas, adicionar o confundimento que deveriam remover
  ])
  #v(0.2em)
  #block(text(size: 11.5pt)[Ruído de estimação, model shift, e o alcance de uma penalidade])
  #v(1.1em)
  #text(size: 10.5pt)[Luan Ivepe]
  #v(0.2em)
  #text(size: 9.5pt, style: "italic")[Pesquisador independente]
  #v(0.2em)
  #text(size: 9pt)[#link("mailto:luanivepe@gmail.com")[luanivepe\@gmail.com]]
  #v(0.2em)
  #text(size: 9pt)[Rascunho de preprint, revisado após parecer — #datetime.today().display("[day]/[month]/[year]")]
  #v(0.4em)
  #text(size: 8.5pt, style: "italic")[Versão em português. A versão de referência, para submissão, é a inglesa.]
]

#v(1.2em)

#block(inset: (x: 1.2em), [
  #text(weight: "bold")[Resumo] #h(0.6em)
  A aceleração de idade epigenética em sangue é rotineiramente corrigida pela
  composição de células imunes, regredindo a idade do relógio sobre as proporções
  celulares estimadas. Dentro da coorte onde é ajustada, a correção faz o que
  promete; estudamos o que acontece quando seus coeficientes são aplicados a outra
  coorte. Em seis coortes públicas de sangue, pontuando só relógios que nunca
  treinaram nas coortes envolvidas, uma correção ajustada em quarenta amostras e
  transportada aumentou o sinal de composição em 83% a 95% dos sorteios nos
  relógios que estimam idade, conforme a coorte de ajuste; no relógio de ritmo de
  envelhecimento DunedinPACE foi neutra. O erro tem duas partes. O ruído de
  estimação cai com $1\/n$ e é reproduzido quase exatamente por coeficientes sem
  informação. O model shift — um efeito de composição que difere entre coortes —
  não caiu com o tamanho de ajuste até 2.639 amostras e responde pelos piores
  transportes. Como depende dos coeficientes do próprio alvo, nada calculável de
  antemão certificou um transporte como seguro: um previsor em forma fechada
  chamou de seguros 30 de 73 transportes nocivos. Uma penalidade ridge fixa
  reduziu os transportes nocivos (par × relógio) de 23 de 72 para 1; uma
  penalidade que diminui com o tamanho amostral, ou escolhida por validação
  cruzada, não. Na saliva, onde a composição responde por até metade da
  aceleração, os transportes entre coortes foram nocivos em 7 de 8 células com ou
  sem penalidade (até +152%), embora a correção funcionasse dentro de cada coorte,
  e reunir estudos não ajudou. A diferença é o sinal: o encolhimento ajuda quando
  o coeficiente do alvo tem o mesmo sinal, o que se manteve entre as coortes de
  sangue e falhou na saliva, onde a inclinação de um relógio na fração imune se
  inverteu entre coortes.
])

= Introdução

A composição do sangue muda com a idade: linfócitos diminuem, células mieloides
aumentam, e células de memória substituem as naive. Cada subtipo de leucócito tem
seu perfil de metilação, então um relógio aplicado a sangue total mede em parte a
composição da amostra @jaffe2014. Com doze tipos celulares resolvidos, a
composição explica de 13% a 34% da variância de aceleração de idade, conforme o
relógio @zhang2024, e células T CD8 naive leem 15 a 20 anos mais novas que CD8
de memória efetora do mesmo doador @tomusiak2024.

O remédio usual é regredir a idade do relógio sobre a idade cronológica e as
proporções estimadas, e guardar o resíduo. Isso costuma ser feito dentro do
próprio conjunto de dados, onde remove exatamente o termo linear de composição.
Pelo mesmo motivo não pode ser verificado ali: um resíduo de mínimos quadrados é
ortogonal aos preditores qualquer que seja a qualidade dos coeficientes. As
avaliações publicadas de ajuste por tipo celular são por simulação e dentro do
conjunto de dados @mcgregor2016.

Examinamos o caso em que os coeficientes saem da coorte que os produziu:
reutilização de coeficientes publicados, coortes pequenas que os tomam de
coortes maiores, ou uma correção fixa aplicada a amostras novas. Testamos isso no sangue e na
saliva, onde a composição pesa mais e uma adaptação publicada já a transporta.

= Métodos

== Coortes, painéis de referência e relógios

Quatro séries públicas de sangue total em 450k com idade cronológica: GSE40279
($n = 656$, idades 19–101), GSE61151 ($n = 184$), GSE50660 ($n = 464$, coorte de
tabagismo) e GSE42861 ($n = 689$, caso-controle de artrite reumatoide). As
proporções foram estimadas por mínimos quadrados não negativos com restrição,
contra um painel de seis tipos (GSE35069) e um de doze tipos que separa
linfócitos naive e de memória (GSE167998) @salas2022, ambos construídos aqui; o
painel de doze recupera proporções conhecidas de misturas com $r$ de 0,79 e erro
absoluto médio de 0,027.

Relógios: Horvath 2013 @horvath2013, Levine 2018 @levine2018 e Horvath 2018
@horvath2018. Um relógio é excluído de todo par que envolve uma coorte em que ele
treinou. Hannum 2013 e Horvath 2013 treinaram no GSE40279 (o Horvath 2013 o lista
como conjunto de treino 3; o GSE42861 foi só de teste). Levine 2018 (InCHIANTI) e
Horvath 2018 não treinaram em nenhuma das quatro. Os pares com GSE40279 são
pontuados só com Levine 2018 e Horvath 2018. Uma quinta coorte, o GSE132203
($n = 795$, array EPIC, majoritariamente afro-americana), testa a replicação em
outra geração de array e outra ancestralidade; nela só se pontuam relógios com
cobertura de sondas acima de 95% (Levine 2018, Horvath 2018). Como relógio de outro tipo,
acrescentamos o DunedinPACE @belsky2022, que estima o ritmo de envelhecimento e
foi treinado numa coorte fora do GEO; nós o reimplementamos a partir dos dados de
modelo publicados no pacote e o validamos (médias por coorte de 0,93 a 1,05;
fumantes atuais 0,14 mais rápidos que quem nunca fumou, $p = 3 times 10^(-11)$).
Uma sexta coorte, GSE55763 @lehne2015 (450k, Londres), serve de coorte de ajuste
grande: sem os 72 arrays de réplica técnica, tem 2.639 adultos não aparentados
(24 a 75 anos); é posterior ao Horvath 2013, não está no treino de nenhum relógio
e os quatro relógios de idade passaram nela nas checagens de cobertura e idade.

A saliva, mistura de epitélio bucal e leucócitos, foi testada em três coortes
adultas: GSE232891 (EPIC, 552 pessoas; doença inflamatória intestinal e
controles), GSE232332 (EPIC, 265 após remover réplicas técnicas; câncer de esôfago
e controles) e GSE78874 (450k, 259; betas calculados do sinal bruto). As
proporções vieram das referências do EpiDISH @teschendorff2017 @zheng2018: um
ajuste hierárquico de nove tipos (epitélio, fibroblasto, sete subtipos imunes) e
uma medição de três tipos (epitélio, fibroblasto, imune), que compartilham o
primeiro passo e por isso favorecem a correção na medição. As duas primeiras
coortes vêm do mesmo grupo e seus arquivos não trazem sondas de genotipagem, então
não foi possível excluir pessoas em comum e elas nunca foram pareadas. Levine 2018
e Horvath 2018 passaram nas checagens de cobertura e idade nas três.

== Correção e pontuação

Na coorte de ajuste, ajustamos $y = beta_0 + beta_1 a + C beta_c$ (idade do
relógio $y$, idade cronológica $a$, proporções $C$ com uma coluna descartada) e
aplicamos $y - (C_"teste" - macron(C)_"ajuste") beta_c$ na coorte de teste. O
desfecho é o termo de composição que sobra no resíduo de idade da coorte de
teste — o incremento de $R^2$ da composição sobre a idade cronológica, menos um
nulo de permutação, como fração da variância de aceleração não corrigida. O dano
líquido $Delta$ é essa fração depois da correção menos antes; $Delta$ positivo
significa que a correção transportada foi pior do que nenhuma. A correção é
ajustada com doze tipos e pontuada com seis, para não ser avaliada pela própria
representação; análises de sensibilidade pontuam com doze tipos e só com as
colunas naive/memória.

Toda contagem é por relógio: um par em que um relógio é prejudicado e outro
ajudado conta como duas células. Agregar relógios dentro do par escondeu boa
parte do dano numa versão anterior desta análise.

== Dois componentes do erro de transporte

Seja $b$ ajustado na coorte A e $S_B$ a covariância de composição da coorte B,
ambos depois de retirar intercepto e idade. A composição que sobra em B é
$(beta_B - b)' S_B (beta_B - b)$. Com $b = beta_A + e$ e
$"Cov"(e) = (sigma^2 \/ n) S_A^(-1)$, seu valor esperado é um termo de
especificação $(beta_A - beta_B)' S_B (beta_A - beta_B)$ mais um termo de
estimação $sigma^2 dot tr(S_A^(-1) S_B) \/ n$. Chamamos $tr(S_A^(-1) S_B)\/n$ de
índice de transporte. O termo de estimação é o excesso de risco padrão de mínimos
quadrados sob _covariate shift_ @eyre2024; o de especificação é _model shift_
@lei2021. Não reivindicamos nenhum dos dois como novo.

Uma referência sem informação é obtida embaralhando linhas inteiras da matriz de
composição da coorte de ajuste antes de ajustar, o que preserva a colinearidade
entre tipos. Ela é medida em cada tamanho de ajuste.

== Penalidade, subamostragem e inferência

Os coeficientes ridge são ajustados na composição padronizada e residualizada
pela idade, com penalidade $alpha$ vezes o autovalor médio, para que $alpha$ seja
comparável entre coortes e tamanhos; $alpha = 0$ reproduz mínimos quadrados. Sob
mudança de distribuição a penalidade ótima nem precisa ser positiva
@patil2024; uma penalidade positiva fixa é usada aqui como padrão conservador.
As subamostras são estratificadas por decil de idade. Onde as configurações
compartilham coortes, a significância é avaliada por permutação em blocos
@winkler2015, trocando perfis inteiros de índice entre pares de coortes.

= Resultados

== Uma correção transportada pode ser pior do que nenhuma

#figure(
  image("figures/pt/fig1_curve.png", width: 95%),
  caption: [*Dano líquido de uma correção transportada, por tamanho de ajuste.*
  Ajustada em subamostras estratificadas por idade do GSE40279 e aplicada a três
  coortes externas; Levine 2018 e Horvath 2018; 30 sorteios por tamanho. Acima de
  zero a correção foi pior do que nenhuma. Laranja: o mesmo procedimento com as
  linhas de composição embaralhadas antes do ajuste.],
) <fig1>

Ajustada em 40 amostras e transportada, a correção teve $Delta$ mediano de +18,9%
(IQR de +8,2% a +34,2%) e foi nociva em 93% dos sorteios (@fig1). Em três
conjuntos independentes de sorteios, a mediana em 40 amostras variou de +11,8% a
+18,9%. A mediana ficou perto de zero entre 160 e 240 amostras e em −1,6% nas 656
completas, onde 33% dos valores (relógio × coorte) ainda eram nocivos. Pontuado
com doze tipos em vez de seis, o dano em 40 amostras foi de +48,1%; só nas
colunas naive/memória, +30,6%. Essa medição usa o mesmo painel do ajuste e é
enviesada a favor da correção, então a curva de seis tipos subestima o dano.

== Dois componentes

#figure(
  image("figures/pt/fig2_components.png", width: 95%),
  caption: [*Composição deixada pelos coeficientes reais e pelos embaralhados.* O
  ajuste real mantém um piso de cerca de 4,5 pontos em tamanho cheio; acima desse
  piso, o excesso acompanha a referência embaralhada.],
) <fig2>

Coeficientes embaralhados não carregam informação, e mesmo assim deixaram +15,6%
de dano líquido em 40 amostras, caindo com o tamanho a uma inclinação log-log de
−0,93 (@fig2). A maior parte do dano em amostras pequenas é, portanto, ruído de
estimação. O ajuste real, porém, ainda deixou 4,5 pontos de composição em 656
amostras, onde a referência embaralhada deixou 1,0. Descontado esse piso, o resto
acompanhou de perto a referência embaralhada (9,6 contra 9,6 pontos em 60
amostras; 3,1 contra 3,1 em 120; 1,2 contra 1,3 em 240).

O piso não é artefato de pontuar com outro painel: aplicada dentro do GSE40279, a
mesma correção removeu 95% (Levine) e 85% (Horvath 2018) do sinal de composição
de seis tipos. Transportada em tamanho cheio, removeu de 92% a −109% conforme o
par — no pior caso, dobrando o sinal.

== Model shift

#figure(
  image("figures/pt/fig3_model_shift.png", width: 95%),
  caption: [*Termo de especificação contra a composição que sobra em tamanho
  cheio.* 24 configurações (par direcionado × relógio). O termo é corrigido pelo
  ruído de estimação dos dois ajustes.],
) <fig3>

O termo de especificação corrigido ordenou a sobra em tamanho cheio com $rho$ de
Spearman de 0,633 ($p$ = 0,0009; @fig3). Ele explica o pior transporte
encontrado, do GSE61151 para a coorte de artrite: previsto +26,6% e +51,8% para
os dois relógios, observado +30,8% e +30,5%. Os testes de Wald par a par de
coeficientes iguais rejeitaram em 4 de 12 comparações após Bonferroni, abaixo da
barra pré-definida de 6 (9 de 12 com $p$ nominal abaixo de 0,05). As duas medidas
pesam as diferenças de coeficiente de modos distintos: o Wald pela precisão da
estimativa, o termo de especificação pela variância na coorte-alvo. Dentro das
coortes, os coeficientes do Levine 2018 diferiram entre casos de artrite e
controles e entre quem já fumou e quem nunca fumou ($p$ = 0,013 cada); para o
Horvath 2018 esse teste não rejeitou ($p$ = 0,61 e 0,48). Ajustar por doença e
tabagismo deixou as diferenças entre coortes praticamente inalteradas.

O model shift não é efeito de lote. Dentro do GSE42861, onde casos e controles
compartilham estudo, array e laboratório, ajustamos a correção numa metade
aleatória dos controles e a aplicamos tanto à outra metade quanto aos casos de
artrite (30 divisões pareadas). Aplicada aos casos, ela deixou 7,4, 7,4 e 2,6
pontos a mais de composição do que aplicada a controles (Horvath 2013, Levine
2018, Horvath 2018), embora o índice de transporte fosse menor para os casos —
então nada do excesso é ruído de estimação, e o componente de model shift
estimado foi de 8,7, 8,4 e 6,3 pontos. O Horvath 2018 é afetado mesmo sem o teste
de Wald acima detectar diferença de coeficientes. A escolha da população de
referência também muda o efeito estimado da doença: no Horvath 2018, o efeito da
artrite ajustado por idade foi de −0,74 ano com a correção na coorte inteira e de
−1,39 ano com a correção ajustada só nos controles.

== O que dá para saber antes de transportar

#figure(
  image("figures/pt/fig4_index.png", width: 95%),
  caption: [*Índice de transporte contra dano líquido, por relógio.* Laranja:
  transportes que o previsor em forma fechada rotulou como seguros e que foram
  nocivos. A forma do marcador indica a coorte de ajuste.],
) <fig4>

Como o índice escala com $1\/n$, parte da sua correlação com o dano vem só do
tamanho amostral. Nos três tamanhos disponíveis para todos os pares, uma
permutação em blocos que preserva os tamanhos de cada par deu nulo com mediana
de $rho$ igual a 0,42; o observado foi 0,579 ($p$ = 0,036). Dentro de um único
tamanho de ajuste, o índice ordenou o dano com $rho$ de 0,41, 0,36 e 0,33. Ele
traz informação real, mas modesta, sobre quais pares são arriscados (@fig4).
Contados por relógio, 12 de 40 transportes com índice abaixo de 0,05 foram
nocivos; nenhum limiar do índice delimita uma região segura.

Uma estimativa em forma fechada do dano líquido que supõe coeficientes
compartilhados, $2 hat(sigma)^2 dot "índice" - b' S_B b$, precisa só da coorte de
ajuste e das proporções do alvo. Ela acertou o sinal em 71% das configurações
($rho$ de 0,633), contra 81% ($rho$ de 0,884) de um oráculo que conhece os
coeficientes do próprio alvo. Ela erra para o lado tranquilizador: 30 dos 73
transportes que rotulou como seguros foram nocivos, concentrados onde o termo de
especificação é grande.

== Uma penalidade remove a maior parte do dano

#figure(
  image("figures/pt/fig5_penalty.png", width: 95%),
  caption: [*Transporte sem penalidade contra transporte com penalidade ridge,
  por célula (par × relógio), em tamanho casado* $n = min(n_A, n_B)$.],
) <fig5>

Sem penalidade, 15 de 30 células (par × relógio) foram nocivas. Com $alpha = 3$,
restou uma, em +0,2% (@fig5); com $alpha = 10$, nenhuma. A penalidade tem custo:
onde a correção sem penalidade já ajudava, o benefício mediano caiu de −4,0% para
−3,3% com $alpha = 3$ e para −1,5% com $alpha = 10$. A penalidade fixa também se
sustentou quando as duas coortes foram deconvoluídas com painéis de referência
diferentes (14 células nocivas sem penalidade, 1 com $alpha = 3$) e quando
pontuada com doze tipos ou no eixo naive/memória (nenhuma célula nociva com
$alpha = 3$).

Escolher a penalidade por validação cruzada leave-one-out na coorte de ajuste
não funcionou tão bem: ela escolheu $alpha$ de 0,3 na célula mediana e deixou 30%
das células nocivas, contra 3% com $alpha = 3$ fixo. A validação cruzada otimiza o
ajuste dentro da coorte de ajuste e não enxerga onde os coeficientes serão usados.

O tamanho da penalidade é definido em relação à variância de composição da coorte
de ajuste, então o encolhimento proporcional não diminui com $n$. Uma penalidade
convencional de $lambda$ fixo diminui; igualada a $alpha = 3$ com 40 amostras, ela
vale 0,045 com 2.639. Nos 30 pares direcionados entre as seis coortes em tamanho
pareado (72 células de relógios de idade), ela deixou 11 células nocivas, contra 1
com $alpha = 3$ (23 sem penalidade); as que escaparam foram os transportes com
model shift. Isso bate com a distinção entre covariate shift e regression shift na
regularização ridge ótima @patil2024: um erro que não diminui com $n$ não é contido
por uma penalidade que diminui.

== Um relógio de ritmo de envelhecimento

No DunedinPACE, a correção transportada ajustada em 40 amostras do GSE40279 foi
neutra (mediana de +0,3%, nociva em 51% dos sorteios), embora os coeficientes
embaralhados ainda tenham causado +4,9% de dano: o componente de ruído estava
presente, mas os coeficientes reais removeram sinal genuíno suficiente para
compensá-lo. Ajustada na quinta coorte, a correção foi neutra de novo (−1,5%,
nociva em 42%), então isso parece propriedade do relógio, e não da coorte de
ajuste. O resto se repetiu: em tamanhos casados, 6 de
12 pares direcionados foram nocivos sem penalidade e nenhum com $alpha = 3$, e
uma correção ajustada em controles deixou 10,5 pontos a mais de composição em
casos de artrite do que em outros controles.

== Replicação em outro array e outra ancestralidade

Ajustada em 40 amostras do GSE132203 (EPIC, majoritariamente afro-americana) e
transportada para as quatro coortes de 450k, a correção foi nociva para os
relógios de idade em 95% dos sorteios (mediana de +24,3%; Levine 2018 91%,
Horvath 2018 99%), com os coeficientes embaralhados em +16,4%. Nos 8 pares
direcionados que envolvem essa coorte, 7 de 24 células (par × relógio) foram
nocivas sem penalidade e nenhuma com $alpha = 3$.

== Uma coorte de ajuste quatro vezes maior

Ajustada no GSE55763 e transportada para as outras cinco coortes (13 células de
relógios de idade), a correção foi nociva com 40 amostras em 83% dos sorteios
(mediana de +7,6%). Entre 656 amostras e as 2.639 completas, a composição que ela
deixou caiu só de 1,9 para 1,3 ponto, razão de 0,72 contra os 0,25 que o ruído
em $1\/n$ prevê; a referência embaralhada caiu de 0,5 para 0,1 (inclinação
log-log de −1,03). O piso, portanto, persiste na mediana. Não é uniforme: no
tamanho cheio foi de +4,6 pontos na coorte de artrite e +3,0 no GSE40279, mas
zero no GSE50660 e no GSE61151, e superou a referência embaralhada em 8 de 13
células, abaixo das 9 que fixamos antes. O model shift é propriedade do par de
coortes. Com tantas amostras de ajuste, a correção sem penalidade ajudou em 10 de
13 células, e a penalidade fixa custou benefício (mediana de −3,0% contra −4,0%)
ao eliminar as três células nocivas. Restringir o GSE40279 à faixa etária da
coorte de ajuste (24 a 75 anos) não mudou seu piso (de +3,5 para +3,3 pontos no
Horvath 2018), então a extrapolação de idade não o explica. O piso se concentra
na metade de ancestralidade europeia da coorte (+4,3 e +3,9 pontos, contra −0,4 e
−0,2 na metade hispânica), que nessa coorte também foi processada em placas
separadas; centrar por placa não o removeu.

== Saliva

#figure(
  image("figures/pt/fig6_saliva.png", width: 95%),
  caption: [*Inclinação de cada relógio na fração imune em quatro coortes de
  saliva.* Idade do relógio regredida na idade cronológica e na fração imune;
  inclinação por 10 pontos percentuais, IC de 95%. GSE149747 na linha de base.],
) <fig6>

Antes da correção, a composição respondia por 10,8% a 51,6% da variância da
aceleração nos relógios de idade na saliva, e por 42% a 56% no DunedinPACE.
Ajustada em 40 amostras e transportada entre coortes de saliva, a correção foi
nociva em 67% dos sorteios (mediana de +20,9%; embaralhada, +6,4%). Em tamanho
pareado (259 amostras), 7 de 8 células (par × relógio) foram nocivas sem
penalidade, 7 de 8 com $alpha = 3$ e 6 de 8 com a penalidade que diminui; no
Horvath 2018 transportado para o GSE78874, a correção sem penalidade deixou +116%
e +152% a mais de sinal de composição do que encontrou. Dentro de cada coorte,
ajustar numa metade aleatória e aplicar na outra foi benéfico em 5 de 6 células
de relógios de idade (de −15% a −48%); a exceção tinha 132 amostras de ajuste e
pouco sinal inicial. A falha está, portanto, no transporte. Todo par de saliva
também cruza array e pré-processamento, então diferenças técnicas e biológicas
não se separam, mas normalizar o GSE78874 por quantis para a distribuição EPIC
deixou 7 de 8 células nocivas, e o mesmo fez um ajuste de três tipos com matriz
de composição bem condicionada (número de condição ≈ 1, contra 872 a 1.317 com
nove tipos). A causa está no eixo dominante (@fig6): ajustado pela idade, o Horvath 2018
mudou +1,0 e +1,9 ano a cada 10 pontos de fração imune nas coortes EPIC e −0,9 no
GSE78874. O encolhimento reduz uma correção de sinal trocado (de +116% para +28%
com $alpha = 3$), mas não a torna útil. Numa quarta coorte de saliva em EPIC, de
outro grupo (GSE149747, 44 adultos na linha de base), a inclinação foi de −0,82
(IC de 95% de −1,71 a 0,07), do lado da coorte de 450k e não das outras coortes
EPIC, o que pesa contra o array como explicação; pela regra que fixamos antes, a
comparação foi inconclusiva.

Reunir estudos também não resgatou o transporte: ajustada nas outras coortes de
saliva, com efeito fixo por estudo, e aplicada à coorte deixada de fora, a
correção foi nociva em 6 de 8 células (5 de 8 com $alpha = 3$), porque a
inclinação reunida toma o sinal da maioria do conjunto.

O sangue difere no sinal. Nas seis coortes de sangue, as
inclinações dos dois relógios em células T CD8 naive mantiveram um só sinal, e em
neutrófilos nenhum intervalo de 95% ficou do lado oposto ao das demais coortes;
na saliva, o Horvath 2018 teve duas coortes claramente positivas e uma claramente
negativa, com $I^2$ entre coortes de 97% contra 43% a 84% no sangue. Com
$alpha = 3$, as células de saliva do Levine 2018, cuja inclinação manteve o sinal,
caíram para uma mediana de +3,2%; as do Horvath 2018 ficaram em +21,2%.

= Discussão

A correção de composição dentro da coorte não pode ser validada dentro dela; isso
decorre dos mínimos quadrados, não dos dados. Transportada, o dano em amostras
pequenas é sobretudo ruído de estimação — coeficientes embaralhados fazem quase o
mesmo estrago —, enquanto o dano em amostras grandes vem de diferenças no efeito
da composição entre coortes. Esse segundo componente é invisível sem os
coeficientes do próprio alvo, e um alvo grande o bastante para estimá-los poderia
simplesmente ser corrigido dentro de si.

Para a prática, isso sugere quatro coisas. Quando a coorte-alvo é grande o
bastante, ajuste a correção dentro dela. No sangue, quando os coeficientes
precisam ser transportados, penalize-os com uma penalidade fixa e substancial, em
vez de uma escolhida por validação cruzada na coorte de origem ou de uma que
diminui com $n$; com milhares de amostras de ajuste isso custa cerca de um ponto
de benefício. Na saliva, e plausivelmente em qualquer tecido em que a composição
domina, não transporte: nenhuma penalidade e nenhum agrupamento de estudos o
tornaram seguro. E trate qualquer diagnóstico prévio, inclusive o índice de
transporte, como ordenação de risco, não como garantia. A mesma cautela vale
dentro de um único estudo: uma correção ajustada em controles e aplicada a
pacientes é um transporte, e aqui ela mudou o efeito estimado da doença em até
duas vezes.

O alcance da penalidade decorre do sinal do efeito da composição. Encolher em
direção a zero aproxima um coeficiente transportado de qualquer alvo cujo
coeficiente próprio tenha o mesmo sinal, e não ajuda um alvo de sinal oposto;
entre as coortes de sangue o sinal se manteve, e na saliva não.

O transporte já é prática publicada fora do sangue: uma adaptação de um relógio
de sangue para saliva ajusta termos de composição em cerca de 960 amostras
reunidas e os aplica a estudos separados @galkin2021. Foi julgada pela acurácia
contra a idade cronológica, que não mostra a composição que sobra na aceleração;
nas nossas coortes de saliva, reunir estudos não evitou o dano. Heterogeneidade
dependente da composição dentro de uma coorte de saliva também já foi descrita
@chan2026.

= Limitações

Seis coortes adultas de sangue total, cinco em 450k e uma em EPIC, duas
definidas por doença ou exposição, e quatro coortes adultas de saliva cujos pares
de transporte todos cruzam array e pré-processamento; outros tecidos e crianças não foram
testados, e a medição na saliva compartilha a referência com o ajuste. Os painéis de
referência foram construídos aqui com seleção de sondas mais simples que as
bibliotecas publicadas. As medianas em tamanhos pequenos variam entre conjuntos
independentes de sorteios (de +11,8% a +18,9% em 40 amostras), e o próprio dano
em amostra pequena depende do relógio e da coorte de ajuste (neutro para o
DunedinPACE a partir do GSE40279). O valor da
penalidade é específico deste painel e destes relógios. A decomposição supõe um
efeito linear da composição. A pontuação de sensibilidade com doze tipos usa o
mesmo painel do ajuste. Nada disso diz respeito a se relógios epigenéticos medem
envelhecimento biológico; trata de uma correção aplicada a eles.

= Disponibilidade de dados e código

Todas as séries são públicas (GSE40279, GSE61151, GSE50660, GSE42861, GSE132203,
GSE55763, GSE232891, GSE232332, GSE78874, GSE149747, GSE35069, GSE167998). O código de análise, o registro etapa a etapa com toda conclusão
derrubada e os scripts das figuras estão em
#link("https://github.com/KTHimiko/clock-lab")[github.com/KTHimiko/clock-lab]
(a ser tornado público antes da submissão).

#bibliography("refs.bib", title: "Referências", style: "nature")
