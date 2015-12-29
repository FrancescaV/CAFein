/* We are starting a fresh integration; clear GSL step and evolve
 objects. */
gsl_odeiv_step_reset(integrator->step);
gsl_odeiv_evolve_reset(integrator->evolve);

/* Enter evolution loop.  NOTE: we *always* take at least one
 step. */
while (1) {
	
	
    /****************************************************************************
	 ******** save the previous timestep
	 ****************************************************************************/
    REAL8 told = t;
	
    /****************************************************************************
	 ******** Call the integrator and check for failure.	 
	 ******** If the integrator is succesfull the time-step is incremented to 
	 ******** t by some about hUsed.
	 ******** 
	 ******** current interval: told-->t
	 ******** hUsed = t-told	 
	 ****************************************************************************/
    status = gsl_odeiv_evolve_apply(integrator->evolve,
									integrator->control, integrator->step, integrator->sys, &t, tend, &h,
									yinit);
	
    /* Check for failure, retry if haven't retried too many times
	 already. */
    if (status != GSL_SUCCESS) {
		if (retries--) {
			/* Retries to spare; reduce h, try again.*/
			h /= 10.0;
			continue;
		} else {
			/* Out of retries, bail with status code. */
			integrator->returncode = status;
			break;
		}
    } else {
		/* Successful step, reset retry counter. */
		retries = integrator->retries;
    }

    /****************************************************************************
	 ******** Now interpolate across the interval told --> t
	 ******** we will compute the interpolated values at equidistant steps 
	 ******** 
	 ******** tintp = told + (t-told)*theta = told + hUsed*theta, 0 <= theta <= √1	
	 ******** 
	 ******** set tintp = told
	 ******** set deltat = 100 (randomly chosen. I will have the interpolated values at 100 equidistant points)
	 
	 ****************************************************************************/

	/* Note we square to get an absolute value, because we may be integrating t in the positive or negative direction */
    while ( (tintp + deltat)*(tintp + deltat) < t*t) {
	
		/*increment the time at which I interpolate by some amount*/
		tintp += deltat;
		
		/*  We have to compute h = (t-told) because the integrator returns a suggested next h, not the actual stepsize taken. */
		REAL8 hUsed = t - told;
	
		REAL8 theta = (tintp - told)/hUsed;
		
		/* These are the interpolating coefficients for y(t + h*theta) = ynew + i1*h*k1 + i5*h*k5 + i6*h*k6 + O(h^4). */
		REAL8 i0 = 1.0 + theta*theta*(3.0-4.0*theta);
		REAL8 i1 = -theta*(theta-1.0);
		REAL8 i6 = -4.0*theta*theta*(theta-1.0);
		REAL8 iend = theta*theta*(4.0*theta - 3.0);
		
		/* Grab the k's from the integrator state. */
		rkf45_state_t *rkfState = integrator->step->state;
		REAL8 *k1 = rkfState->k1;
		REAL8 *k6 = rkfState->k6;
		REAL8 *y0 = rkfState->y0;
		
		for (i = 0; i < dim; i++) {
			ytemp[i] = i0*y0[i] + iend*yinit[i] + hUsed*i1*k1[i] + hUsed*i6*k6[i];
		}
		
		/* Store the interpolated value in the output array. */
		count++;
		if ((status = storeStateInOutput(&output, tintp, ytemp, dim,
										 &outputlen, count)) == XLAL_ENOMEM) goto bail_out;
    }
	
	/*terminate the integration if i am at the end of the interval*/
	
	
    /* Now that we have recorded the last interpolated step that we
	 could, check for termination criteria. */
    if (!integrator->stopontestonly && t >= tend) break;
	
}