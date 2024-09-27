// tabread2d~
// use a 1d table to stor data representing a 2d table.
// this object interpolate between the data to create an audio signal

/* 
This software is copyrighted by Miller Puckette and others.  The following
terms (the "Standard Improved BSD License") apply to all files associated with
the software unless explicitly disclaimed in individual files:

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above  
   copyright notice, this list of conditions and the following 
   disclaimer in the documentation and/or other materials provided
   with the distribution.
3. The name of the author may not be used to endorse or promote
   products derived from this software without specific prior 
   written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,   
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

// Cyrille Henry 09 2024

#include "m_pd.h"
#include <math.h>

/******************** tabread2d~ ***********************/

static t_class *tabread2d_tilde_class;

typedef struct _tabread2d_tilde
{
    t_object x_obj;
    int x_npoints;
    int x_npointsX, x_npointsY;
    t_word *x_vec;
    t_symbol *x_arrayname;
    t_inlet *x_in2;
    t_outlet *x_out;
    t_float f;
} t_tabread2d_tilde;


void tabread2d_tilde_set(t_tabread2d_tilde *x, t_symbol *s)
{
    t_garray *a;
    
    x->x_arrayname = s;
    if (!(a = (t_garray *)pd_findbyclass(x->x_arrayname, garray_class)))
    {
        if (*s->s_name)
            pd_error(x, "tabread2d~: %s: no such array", x->x_arrayname->s_name);
        x->x_vec = 0;
    }
    else if (!garray_getfloatwords(a, &x->x_npoints, &x->x_vec))
    {
        pd_error(x, "%s: bad template for tabread2d~", x->x_arrayname->s_name);
        x->x_vec = 0;
    }
    else { 
    	garray_usedindsp(a);
    	if (x->x_npoints < x->x_npointsX * x->x_npointsY) {
   		pd_error(x, "not enought points for %d*%d interpolation", x->x_npointsX, x->x_npointsY);
   		x->x_vec = 0;
   		}
   	}	
}

void tabread2d_tilde_size(t_tabread2d_tilde *x, t_float sizex, t_float sizey)
{
    x->x_npointsX = sizex;
    x->x_npointsY = sizey;
}

void *tabread2d_tilde_new(t_symbol *s, t_floatarg sizex, t_floatarg sizey)
{
    t_tabread2d_tilde *x = (t_tabread2d_tilde *)pd_new(tabread2d_tilde_class);
    x->x_arrayname = s;
    x->x_vec = 0;
    x->x_out = outlet_new(&x->x_obj, gensym("signal"));
    x->x_in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);

    tabread2d_tilde_set(x, x->x_arrayname);
    tabread2d_tilde_size(x, sizex, sizey);
    return (void *)x;
}

float tabread2d_read (t_tabread2d_tilde *x, int X, int Y) {
	int index = X + Y * x->x_npointsX;
	return x->x_vec[index].w_float;
}

float CubicHermite (float A, float B, float C, float D, float t)
{
    float a = -A / 2.0f + (3.0f*B) / 2.0f - (3.0f*C) / 2.0f + D / 2.0f;
    float b = A - (5.0f*B) / 2.0f + 2.0f*C - D / 2.0f;
    float c = -A / 2.0f + C / 2.0f;
    float d = B;
 
    return a*t*t*t + b*t*t + c*t + d;
}
 
 
t_int *tabread2d_tilde_perform(t_int *w)
{
    t_tabread2d_tilde *x = (t_tabread2d_tilde *)(w[1]);
    t_sample *in1 = (t_sample *)(w[2]);
    t_sample *in2 = (t_sample *)(w[3]);
    t_sample *out = (t_sample *)(w[4]);
    int n = (int)(w[5]);    
    t_word *buf = x->x_vec;
    int i;

    if (!buf) {
    	while (n--) *out++ = 0;
    	return (w+6);
    }
    
    for (i = 0; i < n; i++)
    {
        t_sample findexX = *in1++;
        t_sample findexY = *in2++;
        int iindexX, iindexY;
        float fractX, fractY;
        
		iindexX = (int)findexX;
		fractX = findexX - (int)findexX;
		iindexX =  iindexX % x->x_npointsX;
		if (fractX < 0) { fractX += 1.; iindexX-=1;}
		if (iindexX < 0) iindexX += x->x_npointsX;
		iindexX += x->x_npointsX;// pour etre sur qu'on reste positif malgré le -1

		iindexY = (int)findexY;
		fractY = findexY - (int)findexY;
		iindexY =  iindexY % x->x_npointsY;
		if (fractY < 0) { fractY += 1.; iindexY-=1;}
		if (iindexY < 0) iindexY += x->x_npointsY;
		iindexY += x->x_npointsY;

	 
		float C00 = tabread2d_read(x, (iindexX-1) % x->x_npointsX	, (iindexY-1)%x->x_npointsY );
		float C10 = tabread2d_read(x, iindexX % x->x_npointsX		, (iindexY-1)%x->x_npointsY );
		float C20 = tabread2d_read(x, (iindexX+1)%x->x_npointsX		, (iindexY-1)%x->x_npointsY );
		float C30 = tabread2d_read(x, (iindexX+2)%x->x_npointsX		, (iindexY-1)%x->x_npointsY );
		float C01 = tabread2d_read(x, (iindexX-1) % x->x_npointsX	, iindexY %x->x_npointsY 	);
		float C11 = tabread2d_read(x, iindexX % x->x_npointsX		, iindexY %x->x_npointsY 	);
		float C21 = tabread2d_read(x, (iindexX+1)%x->x_npointsX		, iindexY %x->x_npointsY 	);
		float C31 = tabread2d_read(x, (iindexX+2)%x->x_npointsX		, iindexY %x->x_npointsY 	);
		float C02 = tabread2d_read(x, (iindexX-1) % x->x_npointsX	, (iindexY+1)%x->x_npointsY );
		float C12 = tabread2d_read(x, iindexX % x->x_npointsX		, (iindexY+1)%x->x_npointsY	);
		float C22 = tabread2d_read(x, (iindexX+1)%x->x_npointsX		, (iindexY+1)%x->x_npointsY );
		float C32 = tabread2d_read(x, (iindexX+2)%x->x_npointsX		, (iindexY+1)%x->x_npointsY );
		float C03 = tabread2d_read(x, (iindexX-1) % x->x_npointsX	, (iindexY+2)%x->x_npointsY	);
		float C13 = tabread2d_read(x, iindexX % x->x_npointsX		, (iindexY+2)%x->x_npointsY	);
		float C23 = tabread2d_read(x, (iindexX+1)%x->x_npointsX		, (iindexY+2)%x->x_npointsY	);
		float C33 = tabread2d_read(x, (iindexX+2)%x->x_npointsX		, (iindexY+2)%x->x_npointsY	);
		
		float col0 = CubicHermite(C00, C10, C20, C30, fractX);
        float col1 = CubicHermite(C01, C11, C21, C31, fractX);
        float col2 = CubicHermite(C02, C12, C22, C32, fractX);
        float col3 = CubicHermite(C03, C13, C23, C33, fractX);
        
        *out++ = CubicHermite(col0, col1, col2, col3, fractY);
    }
    return (w+6);

}

void tabread2d_tilde_dsp(t_tabread2d_tilde *x, t_signal **sp)
{
    dsp_add(tabread2d_tilde_perform, 5, x,
        sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
}

void tabread2d_tilde_free(t_tabread2d_tilde *x)
{
  inlet_free(x->x_in2);
  outlet_free(x->x_out);
}

void tabread2d_tilde_setup(void)
{
    tabread2d_tilde_class = class_new(gensym("tabread2d~"),
        (t_newmethod)tabread2d_tilde_new, 
        (t_method)tabread2d_tilde_free,
        sizeof(t_tabread2d_tilde), 
        CLASS_DEFAULT, 
        A_DEFSYM, A_DEFFLOAT, A_DEFFLOAT, 0);
        
    CLASS_MAINSIGNALIN(tabread2d_tilde_class, t_tabread2d_tilde, f);
    class_addmethod(tabread2d_tilde_class, (t_method)tabread2d_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(tabread2d_tilde_class, (t_method)tabread2d_tilde_set,
        gensym("set"), A_SYMBOL, 0);
    class_addmethod(tabread2d_tilde_class, (t_method)tabread2d_tilde_size,
        gensym("size"), A_FLOAT, A_FLOAT, 0);
        
}
