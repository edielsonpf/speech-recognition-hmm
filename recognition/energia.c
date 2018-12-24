//---------------------------------------------------------------------------
// Calcula a energia de um frame
double calcEnergia
(
  double *xw,  // sinal a ser analisado
  int Lwindow // numero de amostras do sinal
)
{
  int register i; // contador
  double soma=0.0; // armazena a energia do sinal

  for (i=0;i<Lwindow;i++)
    soma += xw[i]*xw[i];

  soma /= (double)Lwindow;

  return soma;
}

//----------------------------------------------------------------------------
// Funcao que calcula parametros delta e delta-delta p/ os parametros log energia normalizada
void CalculaDeltaEnergia
(
  double **energia, // vetor com os parametros log-energia, delta e delta-delta
  int n_quadros, // numero de quadros na locucao
  int delta // janelas a esquerda e a direita a serem consideradas para o calculo dos parametros delta
)
{
  int register i,j; // contadores

  // Zerando ponteiro energia
  for (i=0;i<n_quadros;i++)
    for (j=1;j<3;j++)
      energia[i][j] = 0;

  // Calculando parametros delta
  for (i=0;i<n_quadros;i++)
    for (j=-delta;j<=delta;j++)
      if (((i+j)>=0) && ((i+j)<n_quadros))
        energia[i][1] += -(double)j*energia[i+j][0];

  for (i=0;i<n_quadros;i++)
    energia[i][1] /= (2.0*(double)delta + 1.0);

  // Calculando parametros delta-delta
  for (i=0;i<n_quadros;i++)
    for (j=-delta;j<=delta;j++)
      if (((i+j)>=0) && ((i+j)<n_quadros))
        energia[i][2] += -(double)j*energia[i+j][1];

  for (i=0;i<n_quadros;i++)
    energia[i][2] /= (2.0*(double)delta+ 1.0);
}
