/* fft16.c  complex data stored  re, im, re, im, re, ... im */

void fft16(float Z[32]) /* input data points and output  [0] to [31] */
{
static float W[32];    /* scratch vector, used many times, */
static float E[18] =   /* constants for FFT algorithm */
  { 1.0,       0.0,       0.92388,  0.382683,
    0.707107,  0.707107,  0.382683, 0.92388,
    0.0,       1.0,      -0.382683, 0.92388,
    -0.707107,  0.707107, -0.92388,  0.382683,
    -1.0,       0.0};
    float Tre, Tim;
    int i, j, k, l, m;

    m = 8;
    
    for (l=1; l<=m; l=l+l)
    { 
      k = 0;
      i = 0;
      for (j=l; j<=m; j=k+l)
      { 
        for (; i<j; i++)
        {
          /* W[i+k] = Z[i] + Z[m+i]; complex */
          W[2*(i+k)]   = Z[2*i]   + Z[2*(m+i)];
          W[2*(i+k)+1] = Z[2*i+1] + Z[2*(m+i)+1];

          /* W[i+j] = E[k] * (Z[i] - Z[m+i]); complex */
          Tre = Z[2*i]   - Z[2*(m+i)];
          Tim = Z[2*i+1] - Z[2*(m+i)+1];

          W[2*(i+j)]   = E[2*k] * Tre - E[2*k+1] * Tim;
          W[2*(i+j)+1] = E[2*k] * Tim + E[2*k+1] * Tre; 
        }
        k = j;
      }

      l *= 2;

      /* work back other way without copying */
      k = 0;
      i = 0;
      for (j=l; j<=m; j=k+l)
      {
        for (; i<j; i++)
        {
          /* Z[i+k] = W[i] + W[m+i]; complex */
          Z[2*(i+k)]   = W[2*i]   + W[2*(m+i)];
          Z[2*(i+k)+1] = W[2*i+1] + W[2*(m+i)+1];

          /* Z[i+j] = E[k] * (W[i] - W[m+i]); complex */
          Tre = W[2*i]   - W[2*(m+i)];
          Tim = W[2*i+1] - W[2*(m+i)+1];

          Z[2*(i+j)]   = E[2*k] * Tre - E[2*k+1] * Tim;
          Z[2*(i+j)+1] = E[2*k] * Tim + E[2*k+1] * Tre;
        }
        k = j;
      }
    }

    // result is in Z
} /* end fft16 */