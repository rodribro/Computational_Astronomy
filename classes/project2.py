from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import convolve, savgol_filter
from scipy.optimize import curve_fit
from itertools import combinations
import matplotlib.cm as cm

# For synthetic spectra (like M_p5500g4.0z0.00t1.0_a0.00c0.00n0.00o0.00r0.00s0.00_VIS.spec.FITS)

def get_synth_intensity_wavelength(filename):
    hdul_image_synth = fits.open(filename)

    synth_data_1 = hdul_image_synth[1].data
    wavelength_synth = synth_data_1['wavelength']
    flux_synth = synth_data_1['flux']
    normalized_flux_synth = synth_data_1['normalized flux']

    return wavelength_synth, flux_synth, normalized_flux_synth

def get_intensity_wavelength_star(filename,lambda_min,lambda_max):
    '''Dado o ficheiro fits e os limites de comprimento de onda, devolve o array de fluxo e o de comprimento de onda'''
    
    header=fits.getheader(filename)
    
    # Dados do ficheiro fits
    data=fits.getdata(filename)
    ln=len(data)
    
    # Parâmetros
    ref_px=header['CRPIX1']
    val=header['CRVAL1']
    delt=header['CDELT1']
    
    #criar array dos lambdas
    wavelength=np.arange(val+(1-ref_px)*delt,val+(ln-ref_px)*delt,delt)
    
    # Corrigir erro de elemento a menos
    if len(data) > len(wavelength):
        data = data[:-1] 
    elif len(wavelength) > len(data):
        wavelength = wavelength[:-1]  

    #limitar os lambdas aos limites impostos
    mask = (wavelength >= lambda_min) & (wavelength <= lambda_max)
    wv_filtered=wavelength[mask]
    
    #limitar o fluxo aos limites de lambda impostos
    data_filtered=data[mask]
    
    return wv_filtered, data_filtered

# Plot to superimpose spectra
def plot_spectrum(lambda_min, lambda_max, flux_1, wavelength, flux_1_name, x_axis_label=None, y_axis_label=None):

    # Mask to filter values between desired wavelengths
    mask_real = (wavelength >= lambda_min) & (wavelength <= lambda_max)

    #Apply mask to both spectra and wavelength arrays
    wavelength_filtered = wavelength[mask_real]
    flux_1_filtered = flux_1[mask_real]
  

    # Normalize the flux values to the range [0, 1]
    normalized_flux_1 = (flux_1_filtered - flux_1_filtered.min()) / (flux_1_filtered.max() - flux_1_filtered.min())

    # Plot the normalized data
    plt.figure(figsize=(10, 6))
    plt.plot(wavelength_filtered, normalized_flux_1, label=f"{flux_1_name}")
    
    if(x_axis_label == None):
        plt.xlabel("Wavelength (Å)")
        plt.ylabel("Relative Intensity")
    else:
        plt.xlabel(f"{x_axis_label}")
        plt.ylabel(f"{y_axis_label}")  

    plt.legend()
    plt.show()

def fwhm(wavelength, resolvent_power):

    wavelength_avg = np.mean(wavelength)

    return wavelength_avg/resolvent_power


def sigma_gaussian(fwhm):
    const = 2 * np.sqrt(2*np.log(2))

    return fwhm/const

# Perfil instrumental
def P(wavelength, resolvent_power):
    # wavelentgth mean
    wv_mean = wavelength.mean()

    # fwhm
    fwhm_v = fwhm(wv_mean, resolvent_power)

    # sigma
    sigma = sigma_gaussian(fwhm_v)
    
    # Construir a gaussiana
    pixel_size = wavelength[1] - wavelength[0]  # Supor um tamanho de pixel em Å (pequeno para alta resolução)
    xx = np.arange(-6 * sigma, 6 * sigma, pixel_size)  # Intervalo centrado em zero
    gaussian_profile = np.exp(-0.5 * (xx / sigma) ** 2)
    gaussian_profile /= gaussian_profile.sum()  # Normalizar a gaussiana

    return xx, gaussian_profile

# Plot perfil instrumental
def plot_gaussian_P(wavelength, resolvent_power):

    xx, gaussian_profile = P(wavelength, resolvent_power)

    fwhm_value = fwhm(wavelength,resolvent_power)
    sigma = sigma_gaussian(fwhm_value)
    # Calculate FWHM
    FWHM = 2 * np.sqrt(2 * np.log(2)) * sigma
    half_max_intensity = np.max(gaussian_profile) / 2

    # Plotar a gaussiana
    plt.figure(figsize=(8, 5))
    plt.plot(xx, gaussian_profile, label=f"Gaussiana (σ = {sigma:.2f} Å)", color="blue")
    plt.axvline(0, color="red", linestyle="--", label="Centro (λ médio)")
    plt.axhline(half_max_intensity, color="green", linestyle="--", label=f"FWHM = {FWHM:.2f} Å")
    plt.xlabel("Wavelength (A)")
    plt.ylabel("Intensidade Normalizada")
    plt.title("Perfil Gaussiano (Instrumental)")
    plt.legend()
    plt.grid(True)
    plt.show()

# Função de rotação
def G(wavelength, epsilon, vsinI):
    """
    Calculate the rotational broadening profile G(Δλ).
    
    Parameters:
        vsinI: float
            Projected rotational velocity (in the same units as c).
        wavelength: array-like
            Wavelength array.
        epsilon: float
            Linear limb-darkening coefficient (0 to 1).
    
    Returns:
        array-like: Rotational profile G(Δλ).
    """
    # Vel.luz
    c=299792458.0 #m/s

    #Fatores eq
    central_wavelength = np.mean([wavelength.min(), wavelength.max()])
    delta_lambda_m = (central_wavelength * vsinI) / c

    # Step para intervalo
    pixel_size = wavelength[1] - wavelength[0]

    # Intervalo
    xx=np.arange(-delta_lambda_m,delta_lambda_m, pixel_size)

    # Fatores formula
    fator = (xx) / delta_lambda_m
    
    # Termos 
    term1 = 2 * (1 - epsilon) * np.sqrt(1 - fator**2)
    term2 = (np.pi * epsilon / 2) * (1 - fator**2)
    denominator = np.pi * delta_lambda_m * (1 - epsilon / 3)

    # Final result
    G_value = (term1 + term2) / denominator
    G_value = G_value/G_value.sum() #normalize
    
    return xx, G_value

# Plot perfil rotação
def plot_rotational_profile(xx, g_values):

    plt.title("Rotational Broadening Profile")
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Rotational Profile G(Δλ)")
    plt.plot(xx,g_values)
    plt.show()


# Plot to superimpose spectra
def plot_spectra(lambda_min, lambda_max, flux_1, flux_2, wavelength, flux_1_name, flux_2_name, x_axis_label=None, y_axis_label=None):

    # Mask to filter values between desired wavelengths
    mask_real = (wavelength >= lambda_min) & (wavelength <= lambda_max)

    #Apply mask to both spectra and wavelength arrays
    wavelength_filtered = wavelength[mask_real]
    flux_1_filtered = flux_1[mask_real]
    flux_2_filtered = flux_2[mask_real]

    # Normalize the flux values to the range [0, 1]
    normalized_flux_1 = (flux_1_filtered - flux_1_filtered.min()) / (flux_1_filtered.max() - flux_1_filtered.min())
    normalized_flux_2 = (flux_2_filtered - flux_2_filtered.min()) / (flux_2_filtered.max() - flux_2_filtered.min())

    # Plot the normalized data
    plt.figure(figsize=(10, 6))
    plt.plot(wavelength_filtered, normalized_flux_1, label=f"{flux_1_name}")
    plt.plot(wavelength_filtered, normalized_flux_2, label=f"{flux_2_name}")
    if(x_axis_label == None):
        plt.xlabel("Wavelength (Å)")
        plt.ylabel("Relative Intensity")
    else:
        plt.xlabel(f"{x_axis_label}")
        plt.ylabel(f"{y_axis_label}")  

    plt.legend()
    plt.show()


# Simulate noisy spectrum
def simulate_noisy_spectrum(wavelength, flux, resolution, vsinI, SNR, epsilon=0.6):
    
    
    # Apply rotational broadening
    xx, G_profile = G(wavelength, epsilon, vsinI)
    flux_rotated = convolve(flux, G_profile, mode='same')
    
    # Apply instrumental broadening
    xx, gaussian_kernel = P(wavelength, resolution)
    flux_broadened = convolve(flux_rotated, gaussian_kernel, mode='same')
    
    # Add noise to achieve the desired SNR
    signal = flux_broadened
    noise_std = np.mean(signal) / SNR
    noise = np.random.normal(0, noise_std, size=signal.shape)
    flux_noisy = signal + noise
    
    return wavelength, flux_noisy, flux_broadened

def linear_function(x, m, b):
    return m * x + b

def growth_curve_og_opt_params(multiplets, dataframe):
    
    opt_params = []
    xx_list = []
    og_x = []
    og_y = []
    energy_pot = []
    central_wavelengths = []

    for i in range(len(multiplets)):


        # Select multiplets and save original x and y
        multiplet = dataframe[(dataframe['mult'] >= multiplets[i][0]) & (dataframe['mult'] <= multiplets[i][1])]
        y = np.log(multiplet['W'] / (multiplet['lambda']/1000))
        x = np.array(multiplet['log_gf'])# + np.log(multiplet['lambda']/1000)
        og_x.append(np.array(x))
        og_y.append(np.array(y))

        # List of possible central wavelengths
        central_wavelengths.append(np.array(multiplet['lambda']))

        # Energy Potentials for a multiplet
        energy_potential = np.array(multiplet['EP'])
        energy_pot.append(energy_potential)
        
        # Save min and max of xx to then plot growth curve
        xx_list.append((min(x), max(x)))

        # Optimal parameters and covariance matrix
        popt, pcov = curve_fit(linear_function, x, y, p0=[5,15])
        m_optimal, b_optimal = popt

        opt_params.append((m_optimal, b_optimal))

    return opt_params, xx_list, og_x, og_y, energy_pot, central_wavelengths


def plot_growth_curves(optimal_params, xx_list, og_x, og_y):
    
    cmap = cm.get_cmap('tab20', len(optimal_params))
    for i in range(len(optimal_params)):
        color = cmap(i)
        m, b = optimal_params[i]
        x_min, x_max = xx_list[i]
        x_model = np.linspace(x_min, x_max, 100)
        y_model = linear_function(x_model, m, b)
        orig_x = og_x[i]
        orig_y = og_y[i]
        plt.plot(orig_x, orig_y, '*', color=color)
        plt.plot(x_model, y_model, color=color)
    plt.title("Growth Curves")
    plt.ylabel(r'log(W$\lambda$/$\lambda$)')
    plt.xlabel(r'log($gf\lambda$)')
    plt.show()


def pairwise_excitment_temperature(optimal_params, energy_potential):

    delta = np.sqrt((optimal_params[0][0] - optimal_params[1][0])**2 + (optimal_params[0][1] - optimal_params[1][1])**2) #Δ= sqrt( (m1​−m2​)^2+(b1​−b2​)^2 )
    ep_1 = np.mean(energy_potential[0])
    ep_2 = np.mean(energy_potential[1])

    t_exc = np.abs(5040 * (ep_1 -ep_2)) / delta
    
    return t_exc


def find_solar_temperature_combinations(multiplets, dataframe, solar_temp=5780, tolerance=150):
    """
    Encontra combinações de multipletos que geram uma temperatura de excitação próxima à temperatura solar.
    
    Parâmetros:
        multiplets (list): Lista de multipletos.
        dataframe (pd.DataFrame): Dados contendo colunas 'mult', 'W', 'lambda', 'log_gf', e 'EP'.
        solar_temp (float): Temperatura solar alvo (K). Default: 5772 K.
        tolerance (float): Diferença máxima aceitável em Kelvin. Default: 100 K.
    
    Retorno:
        valid_combinations (list): Lista de combinações válidas [(multiplet_1, multiplet_2, T_exc), ...].
    """
    valid_combinations = []
    valid_multiplets = set() # Easier to use a set so it only adds unique values
    # Itera sobre todas as combinações de pares de multipletos
    for comb in combinations(multiplets, 2):
        # Extrai os dois multipletos da combinação
        mult_1, mult_2 = comb
        
        # Calcula os parâmetros otimizados e potenciais de excitação para os dois multipletos
        opt_params, _, _, _, energy_pot, _ = growth_curve_og_opt_params([mult_1, mult_2], dataframe)
        
        # Calcula a temperatura de excitação para o par
        t_exc = pairwise_excitment_temperature(opt_params, energy_pot)
        
        # Verifica se a temperatura está dentro do intervalo aceitável
        if abs(t_exc - solar_temp) <= tolerance:
            valid_multiplets.add(mult_1)
            valid_multiplets.add(mult_2)
            valid_combinations.append((mult_1, mult_2))#, t_exc))

    # Transform into list to be able to iterate over it       
    valid_multiplets = list(valid_multiplets)

    return valid_combinations, valid_multiplets


# Definição da função gaussiana
def w_gaussian(lambdas, A, lambda_c, sigma, B):
    return A * np.exp(-((lambdas - lambda_c)**2) / (2 * sigma**2)) + B


def find_nearest_risca(wavelength, intensity, lambda_0, delta_lambda=0.1, threshold=0.05):
    """
    Detecta a risca mais próxima de um comprimento de onda central λ0.
    
    Parameters:
        wavelength (array): Comprimento de onda do espectro.
        intensity (array): Intensidade do espectro.
        lambda_0 (float): Comprimento de onda central estimado da risca.
        delta_lambda (float): Faixa de tolerância ao redor de λ0 (default=0.1 nm).
        threshold (float): Limite mínimo de variação relativa para detectar a risca (default=0.05).
    
    Returns:
        lambda_c (float): Comprimento de onda central da risca detectada.
        risca_found (bool): Indica se a risca foi encontrada.
    """
    # Isolar faixa de interesse ao redor de lambda_0
    mask = (wavelength >= lambda_0 - delta_lambda) & (wavelength <= lambda_0 + delta_lambda)
    wavelength_subset = wavelength[mask]
    intensity_subset = intensity[mask]
    
    # Verificar se a faixa contém dados suficientes
    if len(wavelength_subset) == 0:
        print(f"Nenhuma faixa encontrada ao redor de λ0 = {lambda_0:.2f}")
        return None, False

    # Encontrar mínimo local (risca de absorção)
    min_idx = np.argmin(intensity_subset)
    lambda_c = wavelength_subset[min_idx]
    I_continuo = max(intensity_subset)
    I_min = intensity_subset[min_idx]
    variation = (I_continuo - I_min) / I_continuo

    # Verificar se o mínimo é significativo
    if variation > threshold:
        return lambda_c
    else:
        #print(f"Variação insuficiente para detectar risca em λ0 = {lambda_0:.2f}")
        return None
    
def spectrum_central_wavelengths(central_wavelengths, star_wv, star_flux):
    star_lambda_c = []
    for i in range(len(central_wavelengths)):
        for j in range(len(central_wavelengths[i])):
            risca = find_nearest_risca(star_wv, star_flux, central_wavelengths[i][j])
            star_lambda_c.append(risca)
    
    # Identificar os índices das entradas que são None
    indices_to_remove = [i for i, risca in enumerate(star_lambda_c) if risca is None]

    # Remove None entries that correspond to non-existing central wavelengths
    star_lambda_c = [risca for risca in star_lambda_c if risca is not None]


    return star_lambda_c, indices_to_remove


def fit_gaussian(wavelength, intensity, lambda_c, delta_lambda=0.5):
    """
    Ajusta uma gaussiana a uma risca no espectro.

    Parameters:
        wavelength (array): Comprimento de onda do espectro.
        intensity (array): Intensidade do espectro.
        lambda_c (float): Comprimento de onda central estimado da risca.
        delta_lambda (float): Largura da faixa ao redor de λc para análise.

    Returns:
        params (tuple): Parâmetros ajustados (A, lambda_c, sigma, B).
        W_lambda (float): Largura equivalente da risca.
    """
    # Isolar faixa da risca
    mask = (wavelength >= lambda_c - delta_lambda) & (wavelength <= lambda_c + delta_lambda)
    wavelength_subset = wavelength[mask]
    intensity_subset = intensity[mask]

    # Normalizar o fluxo
    I_continuo = max(intensity_subset)  # Considerar o máximo como contínuo
    intensity_subset_normalized = intensity_subset / I_continuo


    if len(wavelength_subset) == 0:
        raise ValueError(f"Nenhuma faixa encontrada ao redor de λc = {lambda_c:.2f}")

    # Parâmetros iniciais de entrada: [A, lambda_c, sigma, B]
    p0 = [-1, lambda_c, 0.05, 1]
    # Ajuste da gaussiana
    params, _ = curve_fit(w_gaussian, wavelength_subset, intensity_subset_normalized, p0=p0)
    A, lambda_c_fit, sigma, B = params

    # Calcular largura equivalente da risca
    W_lambda = (np.sqrt(2 * np.pi) * np.abs(A) * sigma) / B
    
    return params, W_lambda



def plot_gaussian_fit(wavelength, intensity, lambda_c, delta_lambda=0.5):
    """
    Ajusta uma gaussiana e plota a sobreposição com o espectro da estrela.

    Parameters:
        wavelength (array): Comprimento de onda do espectro.
        intensity (array): Intensidade do espectro.
        lambda_c (float): Comprimento de onda central estimado da risca.
        delta_lambda (float): Largura da faixa ao redor de λc para análise.
    """
    # Isolar a faixa ao redor de λc
    mask = (wavelength >= lambda_c - delta_lambda) & (wavelength <= lambda_c + delta_lambda)
    wavelength_subset = wavelength[mask]
    intensity_subset = intensity[mask]

    # Ajustar a gaussiana
    #p0 = [max(intensity_subset), lambda_c, 0.1, min(intensity_subset)]  # Parâmetros iniciais
    p0 = [-1, lambda_c, 0.05, 1]
    popt, _ = curve_fit(w_gaussian, wavelength_subset, intensity_subset, p0=p0)
    A, mu, sigma, B = popt

    # Gerar a gaussiana ajustada para visualização
    gaussian_fit = w_gaussian(wavelength_subset, *popt)

    # Plotar o espectro original e a gaussiana ajustada
    plt.figure(figsize=(10, 6))
    plt.plot(wavelength_subset, intensity_subset, label='Espectro Original', color='black', alpha=0.7)
    plt.plot(wavelength_subset, gaussian_fit, label='Gaussiana Ajustada', color='red', linestyle='--')
    plt.title('Ajuste da Gaussiana ao Espectro da Estrela')
    plt.xlabel('Comprimento de Onda (Å)')
    plt.ylabel('Intensidade')
    plt.legend()
    plt.show()


def equivalentd_widths(star_lambda_c, wv, flux):
    w_lambda_vals = []
    for i in range(len(star_lambda_c)):
        _ , w_lambda = fit_gaussian(wv, flux, star_lambda_c[i])
        w_lambda_vals.append(w_lambda)
    return w_lambda_vals



def process_multiplet(dataframe, valid_combinations, star_wv, star_flux, index, zero_or_one):

    # Get parameters based on best multiplet combinations
    _, _, org_x, _, energy_potential, central_wv = growth_curve_og_opt_params([valid_combinations[index][zero_or_one]],dataframe)

    # Find star's corresponding central wavelengths
    star_lambda_c, indices_to_rem = spectrum_central_wavelengths(central_wv, star_wv, star_flux)
    
    # Remove non-relevant indices
    if (len(indices_to_rem) != 0):
        org_x = np.delete(org_x, indices_to_rem)
        energy_potential = np.delete(energy_potential, indices_to_rem)

    # Fit gaussian based on the star's central wavelengths
    eq_widths = equivalentd_widths(star_lambda_c, star_wv, star_flux)

    # Calcular log(W_lambda / lambda)
    log_w_by_lambda = np.log10(np.array(eq_widths) / np.array(star_lambda_c))

    # Calcular log(gf * lambda)
    log_gf_by_lambda = np.array(org_x) + np.log(np.array(star_lambda_c)/1000)

    # Ajustar linha reta: y = mx + b
    x = np.array(energy_potential).ravel()
    y = np.array(log_w_by_lambda - log_gf_by_lambda).ravel()

    # Remove nan and inf
    # Create a mask for valid indices (non-NaN and non-Inf)
    valid_mask = ~np.isnan(x) & ~np.isinf(x) & ~np.isnan(y) & ~np.isinf(y)

    # Apply the mask to both x and y to ensure they have the same length
    x = x[valid_mask]
    y = y[valid_mask]

    mean_ep = np.average(x)
    # Debug
    #print(f"[MULT1]\nLarguras Equivalentes 1: {eq_widths}\nLambdas 1: {star_lambda_c}\nlog(gf): {org_x}\nEnergy potentials 1: {energy_potential}\nlog_w_by_lambda: {log_w_by_lambda}\nlog_gf_by_lambda {log_gf_by_lambda}\nx1: {x}\ny1: {y}\n")

    slope, intercept = np.polyfit(x, y, 1)
    #print(f"Slope 1: {slope}\n")



    return slope, intercept, mean_ep

def estimate_excitation_temperature(valid_combinations, dataframe, star_wv, star_flux):
    avg_temps = []


    for i in range(len(valid_combinations)):

        slope_1, intercept_1,  ep_1 = process_multiplet(dataframe,valid_combinations, star_wv, star_flux, i, 0)

        slope_2, intercept_2, ep_2= process_multiplet(dataframe, valid_combinations, star_wv, star_flux, i, 1)

        delta = np.sqrt((slope_1 - slope_2)**2 + (intercept_1 - intercept_2)**2) #Δ= sqrt( (m1​−m2​)^2+(b1​−b2​)^2 )
        t_exc = np.abs(5040 * (ep_1 -ep_2)) / delta

        avg_temps.append(t_exc)
    return avg_temps
