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
  #text(size: 10.5pt)[Luan Ive Pereira]
  #v(0.2em)
  #text(size: 9pt)[#link("https://orcid.org/0009-0000-9343-3742")[ORCID 0009-0000-9343-3742]]
  #v(0.2em)
  #text(size: 9.5pt, style: "italic")[Pesquisador independente]
  #v(0.2em)
  #text(size: 9pt)[#link("mailto:luanivepe@gmail.com")[luanivepe\@gmail.com]]
  #v(0.2em)
  #text(size: 9pt)[Rascunho de preprint — #datetime.today().display("[day]/[month]/[year]")]
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
  relógios que estimam idade; no relógio de ritmo de envelhecimento DunedinPACE
  foi neutra. Nas unidades que um estudo reporta, ela moveu o efeito estimado de
  fumo ou doença em 0,91 ano em relação à estimativa de dentro da coorte, contra
  0,55 ano de não aplicar correção alguma, e inverteu o sinal dele em 3 de 24
  configurações. O erro tem duas partes: o ruído de estimação, que cai com
  $1\/n$ e é reproduzido por coeficientes sem informação, e o model shift — um
  efeito de composição que difere entre coortes —, que não caiu com o tamanho de
  ajuste até 2.639 amostras. Como o model shift depende dos coeficientes do
  próprio alvo, nada calculável de antemão certificou um transporte como seguro:
  de 73 transportes que um previsor em forma fechada chamou de seguros, 30 foram
  nocivos. Uma penalidade ridge fixa reduziu os transportes nocivos (par ×
  relógio) de 23 de 72 para 1, onde uma penalidade que diminui com o tamanho
  amostral, ou escolhida por validação cruzada, não reduziu. Na saliva, onde a
  composição responde por até metade da aceleração, os transportes foram nocivos
  em 7 de 8 células com ou sem penalidade, e reunir estudos não ajudou. O que
  separa os tecidos é quanto as coortes discordam em relação ao efeito médio: o
  desvio-padrão entre coortes da inclinação de composição foi de 0,09 a 0,55 do
  efeito médio no sangue e de 0,91 a 5,34 na saliva. Uma penalidade encolhe o
  coeficiente em direção a zero, o que fica perto do coeficiente de cada coorte só
  no primeiro caso.
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
coortes maiores, ou uma correção fixa aplicada a amostras novas. Coeficientes de ajuste estimados num conjunto de dados e aplicados a outros já são prática publicada: uma adaptação de relógios de sangue para saliva os ajusta em oito estudos reunidos e os aplica a estudos separados @galkin2021, e seu conjunto de treino inclui o GSE78874, uma das coortes em que achamos que o transporte falha.

= Métodos

== Coortes, painéis de referência e relógios

Quatro séries públicas de sangue total em 450k com idade cronológica: GSE40279
($n = 656$, idades 19–101), GSE61151 ($n = 184$), GSE50660 ($n = 464$, coorte de
tabagismo) e GSE42861 ($n = 689$, caso-controle de artrite reumatoide). As
proporções foram estimadas por mínimos quadrados não negativos com restrição,
contra um painel de seis tipos (GSE35069) e um de doze tipos que separa
linfócitos naive e de memória (GSE167998) @salas2022, ambos construídos aqui; o painel de doze recupera proporções conhecidas de
misturas com $r$ de 0,79 e erro absoluto médio de 0,027. Como os dois foram
construídos aqui e compartilham a construção, todo resultado de sangue foi
repetido em bibliotecas publicadas: a referência de doze tipos de Salas et al.
2022 no conjunto de sondas otimizado por IDOL para 450k, no ajuste, e a
referência publicada de sete tipos de sangue, na medição @teschendorff2017, que
compartilham 20 sondas de 600 e 333.

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
(24 a 75 anos); é posterior ao Horvath 2013, não está no treino de nenhum relógio e todos os relógios passaram nela nas checagens de cobertura e idade.

A saliva, mistura de epitélio bucal e leucócitos, foi testada em três coortes
adultas: GSE232891 (EPIC, 552 pessoas; doença inflamatória intestinal e
controles), GSE232332 (EPIC, 265 após remover réplicas técnicas; câncer de esôfago
e controles) e GSE78874 (450k, 259; betas calculados do sinal bruto). As
proporções vieram das referências do EpiDISH @teschendorff2017 @zheng2018: um
ajuste hierárquico de nove tipos (epitélio, fibroblasto, sete subtipos imunes) e
uma medição de três tipos (epitélio, fibroblasto, imune), que compartilham o
primeiro passo e por isso favorecem a correção na medição. As duas primeiras
coortes vêm do mesmo grupo e seus arquivos não trazem sondas de genotipagem, então
não foi possível excluir pessoas em comum e elas nunca foram pareadas. Levine 2018 e Horvath 2018 passaram nas checagens de cobertura e idade nas três. Como essa medição é o
primeiro passo do próprio ajuste, a saliva também foi pontuada com um painel
independente construído a partir do GSE147318 @middleton2022, saliva de crianças
separada em frações imune e epitelial; suas 300 sondas compartilham 3,3% com o
EpiDISH e ele é usado como escore imune relativo, já que sua escala absoluta não
transfere entre estudos. Uma quarta coorte de saliva, GSE149747 (EPIC, de outro grupo), entra com suas 44 amostras de linha de base na comparação de inclinações e nos ajustes reunidos, mas é pequena demais para uma curva de transporte.

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
@patil2024; uma penalidade positiva fixa é usada aqui como padrão conservador. O valor $alpha = 3$ foi fixado nas quatro primeiras coortes, antes de o GSE132203 e o GSE55763 entrarem no projeto, então o comportamento dele nessas duas e na saliva é fora da amostra que o escolheu.
As subamostras são estratificadas por decil de idade. Onde as configurações
compartilham coortes, a significância é avaliada por permutação em blocos
@winkler2015, trocando perfis inteiros de índice entre pares de coortes.

= Resultados

== Uma correção transportada pode ser pior do que nenhuma

#figure(
  image("figures/pt/fig1_curve.png", width: 95%),
  caption: [*Dano líquido de uma correção transportada, por tamanho de ajuste.*
  Ajustada em subamostras estratificadas por idade do GSE40279 e aplicada a três
  coortes externas; Levine 2018 e Horvath 2018; 100 sorteios por tamanho. Acima de
  zero a correção foi pior do que nenhuma. Laranja: o mesmo procedimento com as
  linhas de composição embaralhadas antes do ajuste.],
) <fig1>

Ajustada em 40 amostras e transportada, a correção teve $Delta$ mediano de +16,1 p.p. (IQR de +5,5 a +29,7) e foi nociva em 88% dos sorteios, sobre 100 sorteios por tamanho (@fig1); três conjuntos anteriores de 30 sorteios deram medianas de +11,8 a +18,9, então conjuntos pequenos isolados são instáveis nesse tamanho. A mediana ficou perto de zero entre 160 e 240 amostras e em −1,6 p.p. nas 656 completas, onde 33% dos valores (relógio × coorte) ainda eram nocivos. Pontuado com doze tipos em vez de seis, o dano em 40 amostras foi de +48,1 p.p.; só nas colunas naive/memória, +30,6. Essa medição usa o mesmo painel do ajuste e é
enviesada a favor da correção, então a curva de seis tipos subestima o dano.

== Dois componentes

#figure(
  image("figures/pt/fig2_components.png", width: 95%),
  caption: [*Composição deixada pelo ajuste real, contra um piso mais ruído.* Azul: o que a correção transportada deixa. Cinza: o piso sem ruído, os 3,8 p.p. que sobram em tamanho cheio depois de descontar o que os coeficientes embaralhados ainda fazem ali. Laranja: esse piso mais o dano dos embaralhados em cada tamanho. As duas séries são composição restante, nas mesmas unidades.],
) <fig2>

Coeficientes embaralhados não carregam informação, e mesmo assim fizeram +17,4 p.p. de dano líquido em 40 amostras, caindo com o tamanho a uma inclinação log-log de −1,16 (@fig2). A maior parte do dano em amostras pequenas é, portanto, ruído de estimação. O ajuste real, porém, ainda deixou 4,5 p.p. de composição em 656 amostras, onde os embaralhados ainda faziam 0,8; a diferença, 3,8 p.p., é um piso que o tamanho de ajuste não remove. Somar esse piso ao dano dos embaralhados reproduz a curva do ajuste real dentro de 0,6 p.p. em todos os tamanhos a partir de 60, e a subestima em 1,2 em 40 amostras (@fig2).

O piso não é artefato de pontuar com outro painel: aplicada dentro do GSE40279, a
mesma correção removeu 95% (Levine) e 85% (Horvath 2018) do sinal de composição
de seis tipos. Transportada em tamanho cheio, removeu de 92% a −109% do sinal que encontrou, conforme o par — no pior caso, dobrando-o. (Percentuais do sinal encontrado são relativos; $Delta$ e composição restante estão em p.p. da variância de aceleração ao longo do texto.)

== Model shift

#figure(
  image("figures/pt/fig3_model_shift.png", width: 95%),
  caption: [*Termo de especificação contra a composição que sobra em tamanho
  cheio.* 24 configurações (par direcionado × relógio). O termo é corrigido pelo
  ruído de estimação dos dois ajustes.],
) <fig3>

O termo de especificação corrigido ordenou a sobra em tamanho cheio com $rho$ de
Spearman de 0,633 ($p$ = 0,0009; @fig3). Ele explica o pior transporte em tamanho de ajuste cheio, do GSE61151 para a coorte de artrite: previsto +26,6 e +51,8 p.p. para os dois relógios, observado +30,8 e +30,5. Os testes de Wald par a par de
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
referência também muda o efeito estimado da doença: no Horvath 2018, o efeito da artrite ajustado por idade foi de −0,74 ano (IC de 95% de −1,05 a −0,37) com a correção na coorte inteira e de −1,39 (de −2,08 a −0,80) com a correção ajustada só nos controles, uma diferença de −0,66 (de −1,13 a −0,31) por bootstrap que reajusta as duas correções em cada reamostra.

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
($rho$ de 0,6332 em 126 configurações — por coincidência próximo do $rho$ do
termo de especificação acima, que é outra correlação, em 24), contra 81%
($rho$ de 0,884) de um oráculo que conhece os
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
restou uma, em +0,2 p.p. (@fig5); com $alpha = 10$, nenhuma. A penalidade tem custo:
onde a correção sem penalidade já ajudava, o benefício mediano caiu de −4,0 p.p. para
−3,3 com $alpha = 3$ e para −1,5 com $alpha = 10$. A penalidade fixa também se
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

No DunedinPACE, a correção transportada ajustada em 40 amostras do GSE40279 foi neutra (mediana de +0,3 p.p., nociva em 51% dos sorteios), embora os coeficientes
embaralhados ainda tenham causado +4,9 p.p. de dano: o componente de ruído estava
presente, mas os coeficientes reais removeram sinal genuíno suficiente para
compensá-lo. Ajustada na quinta coorte, a correção foi neutra de novo (−1,5 p.p., nociva em 42%), então isso parece propriedade do relógio, e não da coorte de
ajuste. O resto se repetiu: em tamanhos casados, 6 de
12 pares direcionados foram nocivos sem penalidade e nenhum com $alpha = 3$, e
uma correção ajustada em controles deixou 10,5 pontos a mais de composição em
casos de artrite do que em outros controles.

== Replicação em outro array e outra ancestralidade

Ajustada em 40 amostras do GSE132203 (EPIC, majoritariamente afro-americana) e
transportada para as quatro coortes de 450k, a correção foi nociva para os
relógios de idade em 95% dos sorteios (mediana de +24,3 p.p.; Levine 2018 91%,
Horvath 2018 99%), com os coeficientes embaralhados em +16,4. Nos 8 pares direcionados que envolvem essa coorte, 7 de 24 células (par × relógio) foram nocivas sem penalidade e nenhuma com $alpha = 3$; as 24 contam o DunedinPACE junto dos dois relógios de idade (6 das 16 células de relógios de idade foram nocivas).

== Uma coorte de ajuste quatro vezes maior

Ajustada no GSE55763 e transportada para as outras cinco coortes (13 células de
relógios de idade), a correção foi nociva com 40 amostras em 83% dos sorteios
(mediana de +7,6 p.p.). Entre 656 amostras e as 2.639 completas, a composição que ela
deixou caiu só de 1,9 para 1,3 ponto, razão de 0,72 contra os 0,25 que o ruído
em $1\/n$ prevê; a referência embaralhada caiu de 0,5 para 0,1 (inclinação
log-log de −1,03). O piso, portanto, persiste na mediana. Não é uniforme: no
tamanho cheio foi de +4,6 pontos na coorte de artrite e +3,0 no GSE40279, mas
zero no GSE50660 e no GSE61151, e superou a referência embaralhada em 8 de 13
células, abaixo das 9 que fixamos antes. O model shift é propriedade do par de
coortes. Com tantas amostras de ajuste, a correção sem penalidade ajudou em 10 de
13 células, e a penalidade fixa custou benefício (mediana de −3,0 p.p. contra −4,0)
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
nociva em 67% dos sorteios (mediana de +20,9 p.p.; embaralhada, +6,4). Em tamanho
pareado (259 amostras), 7 de 8 células (par × relógio) foram nocivas sem
penalidade, 7 de 8 com $alpha = 3$ e 6 de 8 com a penalidade que diminui; no Horvath 2018 transportado para o GSE78874, a correção sem penalidade deixou 132 e 168 p.p. de composição onde havia encontrado 16 — de oito a dez vezes o sinal que deveria remover. Dentro de cada coorte,
ajustar numa metade aleatória e aplicar na outra foi benéfico em 5 de 6 células de relógios de idade (de −15 a −48 p.p.); a exceção tinha 132 amostras de ajuste e
pouco sinal inicial. A falha está, portanto, no transporte. Todo par de saliva
também cruza array e pré-processamento, então diferenças técnicas e biológicas
não se separam, mas normalizar o GSE78874 por quantis para a distribuição EPIC
deixou 7 de 8 células nocivas, e o mesmo fez um ajuste de três tipos com matriz
de composição bem condicionada (número de condição ≈ 1, contra 872 a 1.317 com
nove tipos). A causa está no eixo dominante (@fig6): ajustado pela idade, o Horvath 2018 mudou +1,0 e +1,9 ano a cada 10 pontos de fração imune nas coortes EPIC e −0,9 no GSE78874, e o Levine 2018, que mantém o sinal nas quatro coortes, variou de −0,2 a −4,7. O encolhimento leva um coeficiente transportado para zero, então só consegue limitar uma correção que discorda tanto assim do alvo: com $alpha = 3$ a pior célula caiu de +116 para +28 p.p., e as células do Levine 2018, com o sinal intacto, para uma mediana de +3,2 em vez de zero. Numa quarta coorte de saliva em EPIC, de
outro grupo (GSE149747, 44 adultos na linha de base), a inclinação foi de −0,82
(IC de 95% de −1,71 a 0,07), do lado da coorte de 450k e não das outras coortes
EPIC, o que pesa contra o array como explicação; pela regra que fixamos antes, a
comparação foi inconclusiva.

Reunir estudos também não resgatou o transporte: ajustada nas outras coortes de
saliva, com efeito fixo por estudo, e aplicada à coorte deixada de fora, a
correção foi nociva em 6 de 8 células (5 de 8 com $alpha = 3$), porque a
inclinação reunida toma o sinal da maioria do conjunto.

O sangue difere em quanto as coortes discordam em relação ao próprio efeito. Num ajuste de efeitos aleatórios das inclinações por coorte, o desvio-padrão entre coortes dividido pelo efeito médio em módulo foi de 0,09 a 0,55 nas seis coortes de sangue (eixos CD8 naive e neutrófilos, nos dois relógios) e de 0,91 no Levine 2018 e 5,34 no Horvath 2018 na saliva; o $I^2$ foi de 43% a 84% contra 97%. A inversão de sinal é o extremo desse espalhamento, e não é tudo: com o GSE78874 normalizado, a inclinação do Horvath 2018 na saliva é −0,20 (IC de 95% de −0,46 a +0,06), já sem sinal claramente oposto ao das coortes EPIC, e 7 de 8 células seguem nocivas; o Levine 2018 nunca inverte e ainda é nocivo com $alpha = 3$ em 3 de 4 células.

== O que isso faz com uma associação reportada

Composição restante não é o que um estudo reporta. Para quatro exposições — fumo
(GSE50660), artrite reumatoide (GSE42861), doença inflamatória intestinal e
câncer de esôfago (as duas coortes de saliva em EPIC) — reestimamos o coeficiente
da exposição em (idade do relógio) ~ idade + exposição sob três correções:
nenhuma, ajustada dentro do alvo e transportada. Tomando a correção de dentro da
coorte como comparador, uma correção ajustada em 40 amostras de outra coorte
moveu o efeito reportado numa mediana de 0,91 ano, contra 0,55 ano de não aplicar
correção alguma; com $alpha = 3$, 0,48; em tamanho de ajuste cheio, 0,25. De 24
células (alvo × origem × relógio), 13 moveram mais de um ano e 3 inverteram o
sinal; a maior foi de 4,19 anos, no Levine 2018 levado do GSE132203 para a coorte
de artrite, cujo efeito dentro da coorte é de +0,04 ano. Com $alpha = 3$, 6 de 24
ainda moveram mais de um ano.

Essa métrica também mostra uma cauda que a composição restante, limitada por
construção, não mostra: em 40 amostras, 1,2% de 720 sorteios ficaram a mais de 10
anos da estimativa de dentro da coorte, o pior a 37,6 anos. Resolvido como
mínimos quadrados exatos, sem descartar valores singulares próximos de zero,
10,1% passaram de 10 anos, então o tamanho da cauda depende do solucionador,
ainda que a existência dela não dependa; ela se concentra na saliva, onde a
matriz de composição é quase singular. Com $alpha = 3$ nenhum sorteio passou de 10 anos.

A confiabilidade não substitui essa métrica. O GSE55763 mediu 36 pessoas duas
vezes; nesses 72 arrays, que o ajuste nunca viu, o ICC(2,1) técnico da
aceleração caiu depois da correção de dentro do estudo nos quatro relógios de
idade (Horvath 2013 de 0,817 para 0,769), como já se viu que corrigir por
composição reduz a confiabilidade @sehgal2026. Uma correção transportada de 40
amostras o aumentou (0,842), acima do valor sem correção em três relógios. O
motivo é aritmético: a correção quase não mudou a discordância entre os dois
arrays de uma pessoa (de 3,68 para 3,53 no Horvath 2013) e removeu variância
entre pessoas (de 18,32 para 13,56), que é o que uma correção de composição deve
remover. Uma correção que remove menos pontua, portanto, como mais confiável;
a confiabilidade não ordena correções, e o deslocamento de uma associação
reportada ordena.

== Robustez aos painéis de referência

Repetidos nas bibliotecas publicadas, sem nenhuma outra mudança, os resultados de
sangue se mantêm: a correção ajustada em 40 amostras do GSE40279 foi nociva em
77% de 600 sorteios (mediana de +9,1 p.p., contra 88% e +16,1 com os painéis
construídos aqui); nas 60 células (par × relógio) em tamanho pareado, 17 foram
nocivas sem penalidade e nenhuma com $alpha = 3$; e, ajustada em todo o GSE55763,
a composição restante mediana foi de +1,9 p.p. Os dois painéis de doze tipos
concordam no que ambos estimam (neutrófilos com $r$ de 0,985 a 0,997; CD8 naive
de 0,852 a 0,919, nas seis coortes).

Na saliva, o painel independente correlaciona-se de 0,951 a 0,997 com a fração
imune do EpiDISH e vê quase o mesmo sinal de composição antes da correção que a
própria coluna imune do EpiDISH (por exemplo, +1,5 contra −0,1 p.p. no GSE232891
para o Levine 2018; +50,6 contra +50,6 no GSE78874). O que difere é o número de
eixos, não a origem deles: a medição de três tipos vê +18,5 e +31,8 p.p. no
GSE232891, onde o eixo imune sozinho vê cerca de zero, então nas coortes EPIC a
maior parte do sinal de composição na aceleração está nos eixos
epitélio/fibroblasto. Pontuadas só no eixo imune com o painel independente, 6 de
8 células foram nocivas sem penalidade e 5 de 8 com $alpha = 3$, contra 7 e 7 com
a medição compartilhada.

= Discussão

A correção de composição dentro da coorte não pode ser validada dentro dela;
isso decorre dos mínimos quadrados, não dos dados. Ela também não sai de graça:
corrigir por frações imunes dentro de um conjunto de dados reduz a confiabilidade
biológica de quase todo relógio @sehgal2026, então a correção de dentro da coorte
que usamos como comparador é uma referência, não um padrão de verdade. Transportada, o dano em amostras
pequenas é sobretudo ruído de estimação — coeficientes embaralhados fazem quase o
mesmo estrago —, enquanto o dano em amostras grandes vem de diferenças no efeito
da composição entre coortes. Esse segundo componente é invisível sem os
coeficientes do próprio alvo, e um alvo grande o bastante para estimá-los poderia
simplesmente ser corrigido dentro de si.

Nas unidades que um estudo reporta, o custo de pegar emprestado uma correção é de cerca de um ano num efeito de doença ou exposição, e de mais de quatro anos na pior configuração que encontramos — o bastante para mudar o que um artigo conclui.

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

O alcance da penalidade decorre de quanto as coortes discordam em relação ao efeito que se quer corrigir. Encolher em direção a zero troca um coeficiente transportado por um menor, que fica perto do coeficiente próprio de toda coorte quando o espalhamento entre coortes é uma fração do efeito médio, como no sangue, e perto de nenhum quando o espalhamento é tão grande quanto o efeito, como na saliva. A inversão de sinal é o caso extremo, não um mecanismo à parte, e esperaríamos a penalidade falhar sempre que essa razão se aproxima de um.

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
de transporte todos cruzam array e pré-processamento; outros tecidos e crianças
não foram testados. Os painéis construídos aqui foram conferidos contra
bibliotecas publicadas, e a medição na saliva contra um painel independente, sem
mudar nenhuma conclusão. As medianas em tamanhos pequenos são instáveis entre conjuntos independentes de 30 sorteios (de +11,8 a +18,9 p.p. em 40 amostras; a mediana conjunta com 100 sorteios é +16,1), e o próprio dano
em amostra pequena depende do relógio e da coorte de ajuste (neutro para o
DunedinPACE a partir do GSE40279). O valor da
penalidade é específico deste painel e destes relógios. A decomposição supõe um efeito linear da composição; termos quadráticos para os quatro maiores componentes acrescentam uma mediana de 0,002 ao $R^2$ dentro da coorte (significativos em 4 de 12 células coorte × relógio), e uma correção quadrática não transporta melhor que a linear (22 contra 23 células nocivas de 60, mesma mediana). A pontuação de sensibilidade com doze tipos usa o
mesmo painel do ajuste. Nada disso diz respeito a se relógios epigenéticos medem
envelhecimento biológico; trata de uma correção aplicada a eles.

= Disponibilidade de dados e código

Todas as séries são públicas (GSE40279, GSE61151, GSE50660, GSE42861, GSE132203,
GSE55763, GSE232891, GSE232332, GSE78874, GSE149747, GSE147318, GSE35069,
GSE167998); as matrizes de referência publicadas vêm do pacote EpiDISH do
Bioconductor (versão 2.28.0). O código de análise, os scripts das figuras e o
registro etapa a etapa do projeto, com toda conclusão que uma etapa posterior
derrubou, estão em #link("https://github.com/KTHimiko/clock-lab")[github.com/KTHimiko/clock-lab]
(código sob a licença MIT; texto e figuras sob CC BY 4.0). Todo número deste
manuscrito é produzido por um script de lá.

= Declarações

*Uso de IA generativa.* Este estudo foi feito com uso extenso de um sistema de
IA, o Claude (Anthropic), pelo ambiente Claude Code. A IA escreveu o código de
análise, desenhou e rodou as análises com checagens escritas antes de cada
resultado ser lido, buscou e resumiu a literatura e redigiu o texto e as
figuras. O autor definiu os objetivos, decidiu quais perguntas seguir e quais
resultados reportar, obteve leituras críticas informais de versões anteriores e
mandou fazer as correções que elas exigiram, e responde pelo conteúdo. A IA não
é autora.

*Financiamento.* Nenhum.

*Conflito de interesses.* Nenhum.

*Agradecimentos.* Aos leitores das versões anteriores, cujas críticas corrigiram
várias afirmações, cada uma registrada na tabela de correções pública do projeto.

#bibliography("refs.bib", title: "Referências", style: "nature")
