import numpy as np
import matplotlib.pyplot as plt
import kepler_exo as kp
from astropy.time import Time


import matplotlib.transforms as mtransforms

import matplotlib.dates as md
import dateutil

import argparse
import yaml
from scipy.optimize import leastsq


def phase_calc(p, t_obs, period, e, omega):
    t_mod = kp.kepler_Tcent_T0P(period, p, e, omega)
    return t_obs-t_mod


parser = argparse.ArgumentParser(prog='quadrature_check.py', description='quadrature check to easily detect your candidate planets')
parser.add_argument('config_file', type=str, nargs=1, help='config file')
args = parser.parse_args()
file_conf = args.config_file[0]
stream = open(file_conf, 'r')


try:
    config_in = yaml.load(stream, Loader=yaml.FullLoader)
except AttributeError:
    config_in = yaml.load(stream)
    print(' Consider updating YAML')
except:
    print(' Some error happened while reading the configuration file')
    quit()

accepted_extensions = ['.yaml', '.yml', '.conf', '.config', '.input', ]
for extension in accepted_extensions:
    if file_conf.find(extension) > 0:
        yaml_name = file_conf.replace(extension, "")

if 'planets' not in config_in:
    planets = {'b': config_in}
else:
    planets = config_in['planets']

if 'night_start' not in config_in: config_in['night_start'] = '18:00'
if 'night_duration' not in config_in: config_in['night_duration'] = 14.00
if 'quadrature_window' not in config_in: config_in['quadrature_window'] = 0.05
if 'visibility_start' not in config_in: config_in['visibility_start'] = '18:00'
if 'visibility_duration' not in config_in: config_in['visibility_duration'] = 14.00

if len(planets) > 1:
    multiplanets = True
else:
    multiplanets = False

if 'plot_planets' not in config_in: config_in['plot_planets'] = planets


for planet_name, planet in planets.items():

    if 'R' not in planet: planet['R'] = 2.6
    if 'K' not in planet: planet['K'] = 1.000

    if 'quadrature_window' not in planet: planet['quadrature_window'] = config_in['quadrature_window']

    if 'user_phase' not in planet:
        planet['print_phase'] = False
        planet['user_phase'] = [0.6, 0.01]
    else:
        planet['print_phase'] = True

#for planet_name, planet in planets.items():
for planet_name in config_in['plot_planets']:

    planet = planets[planet_name]


    print()
    print(' *** Planet ', planet_name, ' *** ')

    single_orbit_bjd = np.arange(-0.5, 1.5, 0.001, dtype=np.double)
    single_orbit_RV = -1. * np.sin(2*np.pi *  single_orbit_bjd) # Period 1 d, semiapmplitude = 1 m/s

    single_orbit_1period_bjd = np.arange(0.0, 1.0, 0.001, dtype=np.double)
    single_orbit_1period_RV = -1. * np.sin(2*np.pi *  single_orbit_1period_bjd) # Period 1 d, semiapmplitude = 1 m/s

    #upper_quad = [single_orbit_bjd[np.where(single_orbit_RV > 1.0-planet['quadrature_window'])[0][0]],
    #              single_orbit_bjd[np.where(single_orbit_RV > 1.0-planet['quadrature_window'])[0][-1]]]
    phase_upper_quadrature = single_orbit_1period_bjd[np.argmax(single_orbit_1period_RV)]

    upper_quad = [
        single_orbit_bjd[np.where(single_orbit_bjd > phase_upper_quadrature-planet['quadrature_window'])[0][0]],
        single_orbit_bjd[np.where(single_orbit_bjd < phase_upper_quadrature+planet['quadrature_window'])[0][-1]]]

    print('first upper quadrature since reference time --> ', upper_quad)

    #lower_quad = [single_orbit_bjd[np.where(single_orbit_RV < planet['quadrature_window']-1.0)[0][0]],
    #              single_orbit_bjd[np.where(single_orbit_RV < planet['quadrature_window']-1.0)[0][-1]]]
    phase_lower_quadrature = single_orbit_1period_bjd[np.argmin(single_orbit_1period_RV)]
    lower_quad = [
        single_orbit_bjd[np.where(single_orbit_bjd > phase_lower_quadrature-planet['quadrature_window'])[0][0]],
        single_orbit_bjd[np.where(single_orbit_bjd < phase_lower_quadrature+planet['quadrature_window'])[0][-1]]]

    print('first lower quadrature since reference time --> ', lower_quad)

    input_user = [
        single_orbit_bjd[np.where(single_orbit_bjd > planet['user_phase'][0] - planet['user_phase'][1])[0][0]],
        single_orbit_bjd[np.where(single_orbit_bjd < planet['user_phase'][0] + planet['user_phase'][1])[0][-1]]]

    print('first user input phase since reference time --> ', input_user)


    if lower_quad[0] < upper_quad[0]:
        lower_quad = [(i + 1.) for i in lower_quad]

    if input_user[0] < upper_quad[0]:
        order_dict = {0: 'user', 1: 'upper', 2: 'lower'}
        if input_user[1] >= upper_quad[0]:
            input_user[1] = upper_quad[0]
    elif input_user[0] < lower_quad[0]:
        order_dict = {0: 'upper', 1: 'user', 2: 'lower'}
        if input_user[0] <= upper_quad[1]:
            input_user[0] = upper_quad[1]
        if input_user[1] >= lower_quad[0]:
            input_user[1] = lower_quad[0]
    else:
        order_dict = {0: 'upper', 1: 'lower', 2: 'user'}
        if input_user[0] <= lower_quad[1]:
            input_user[0] = lower_quad[1]
    print()

    for date in config_in['dates']:

        date_observations = Time([date+'T'+config_in['night_start']], format='isot', scale='utc')
        observing_interval = [date_observations.jd[0]-0.5, date_observations.jd[0]+config_in['night_duration']/24.+0.5]

        date_visibility = Time([date+'T'+config_in['visibility_start']], format='isot', scale='utc')
        visibility_interval = [date_visibility.jd[0], date_visibility.jd[0]+config_in['visibility_duration']/24.]

        bjd= np.arange(observing_interval[0], observing_interval[1], planet['P']/100.)
        RV = - planet['K'] *  np.sin(2*np.pi* (bjd-planet['Tc'])/planet['P'])


        if multiplanets:
            RV_all = RV * 0.0
            RV_dict = {}
            for pp_name, pp in planets.items():
                RV_all -= pp['K'] * np.sin(2*np.pi* (bjd-pp['Tc'])/pp['P'])
                if planet_name == pp_name:
                    continue
                RV_dict[pp_name] = -pp['K'] * np.sin(2*np.pi* (bjd-pp['Tc'])/pp['P'])

        user_pos = -1

        #fig = plt.figure(figsize=(12, 12))
        fig, ax = plt.subplots(figsize=(8, 8))





        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)

        print(' ---------- night: ', date , ' ---------- ')

        obs_times = open(yaml_name + '_' + date + '_' + planet_name + '.dat', 'w')
        obs_times.write(' ---------- night: '+ date + ' ----------  \n')

        index_ref = [int((i-planet['Tc']) / planet['P']) for i in observing_interval]


        before_rising = [dateutil.parser.parse(date+'T'+config_in['night_start']),
                         dateutil.parser.parse(date+'T'+config_in['visibility_start'])]
        after_setting = [dateutil.parser.parse(date+'T'+config_in['visibility_start']) +
                         dateutil.relativedelta.relativedelta(hours=config_in['visibility_duration']),
                         dateutil.parser.parse(date+'T'+config_in['night_start']) +
                         dateutil.relativedelta.relativedelta(hours=config_in['night_duration'])]

        before_night = [before_rising[0]- dateutil.relativedelta.relativedelta(hours=12), before_rising[0]]
        after_night = [after_setting[1], after_setting[1] + dateutil.relativedelta.relativedelta(hours=12)]
        #date_observations = Time([date+'T'+config_in['night_start']], format='isot', scale='utc')
        #observing_interval = [date_observations.jd[0]-0.5, date_observations.jd[0]+config_in['night_duration']/24.+0.5]

        #date_visibility = Time([date+'T'+config_in['visibility_start']], format='isot', scale='utc')
        #visibility_interval = [date_visibility.jd[0], date_visibility.jd[0]+config_in['visibility_duration']/24.]


        ax.fill_between(before_night, 0, 1, facecolor='black', alpha=0.2, transform=trans, zorder=9)
        ax.fill_between(after_night, 0, 1, facecolor='black', alpha=0.2, transform=trans, zorder=9)
        ax.fill_between(before_rising, 0, 1, facecolor='black', alpha=0.5, transform=trans, zorder=10)
        ax.fill_between(after_setting, 0, 1, facecolor='black', alpha=0.5, transform=trans, zorder=10)



        for ii in range(index_ref[0]-1, index_ref[1]+1):
            upper_quad_day = [(i+ii)*planet['P']+planet['Tc'] for i in upper_quad]
            lower_quad_day = [(i+ii)*planet['P']+planet['Tc'] for i in lower_quad]
            input_user_day = [(i + ii) * planet['P'] + planet['Tc'] for i in input_user]

            t_upper_quad = Time(upper_quad_day, format='jd')
            t_lower_quad = Time(lower_quad_day, format='jd')
            t_input_user = Time(input_user_day, format='jd')

            plot_upper = [dateutil.parser.parse(s) for s in t_upper_quad.iso]
            plot_lower = [dateutil.parser.parse(s) for s in t_lower_quad.iso]
            plot_input = [dateutil.parser.parse(s) for s in t_input_user.iso]

            for i in [0, 1, 2]:
                if order_dict[i] == 'user' and planet['print_phase'] and \
                                input_user_day[1] > visibility_interval[0] and \
                                input_user_day[0] < visibility_interval[1]:
                    obs_times.write('USER INPUT PHASE start: %s ( %f)  end: %s (%f) \n' % (
                        t_input_user.iso[0], input_user_day[0], t_input_user.iso[1], input_user_day[1]))
                    ax.axvline(plot_input[0], c='green')
                    ax.axvline(plot_input[1], c='green')
                    ax.fill_between(plot_input, 0, 1, facecolor='green', alpha=0.5, transform=trans,  label='User-defined phase')

                if order_dict[i] == 'upper' and \
                                upper_quad_day[1] > visibility_interval[0] and \
                                upper_quad_day[0] < visibility_interval[1]:
                    obs_times.write('UPPER QUADRATURE start: %s ( %f)  end: %s (%f) \n' % (
                        t_upper_quad.iso[0], upper_quad_day[0], t_upper_quad.iso[1], upper_quad_day[1]))
                    ax.axvline(plot_upper[0], c='r')
                    ax.axvline(plot_upper[1], c='r')
                    ax.fill_between(plot_upper, 0, 1, facecolor='red', alpha=0.5, transform=trans, label='Upper quadrature')

                if order_dict[i] == 'lower' and \
                                lower_quad_day[1] > visibility_interval[0] and \
                                lower_quad_day[0] < visibility_interval[1]:
                    obs_times.write('LOWER QUADRATURE start: %s ( %f)  end: %s (%f) \n' % (
                        t_lower_quad.iso[0], lower_quad_day[0], t_lower_quad.iso[1], lower_quad_day[1]))
                    ax.axvline(plot_lower[0], c='b')
                    ax.axvline(plot_lower[1], c='b')
                    ax.fill_between(plot_lower, 0, 1, facecolor='blue', alpha=0.5, transform=trans, label='Lower quadrature')

        t_bjd = Time(bjd, format='jd')
        #t_bjd.format = 'isot'

        plot_bjd = [dateutil.parser.parse(s) for s in t_bjd.iso]


        next_day = date_observations.jd[0] + 1.0



        #input_date_ticks = [date+'T18:00', date+'T19:00', date+'T20:00', date+'T21:00', date+'T22:00', date+'T23:00',
        #                    date+'T00:00', date+'T01:00', date+'T02:00', date+'T03:00', date+'T04:00', date+'T05:00']
        #date_ticks = [dateutil.parser.parse(s) for s in input_date_ticks]
        #ax.set_xticks(date_ticks)

        time_start = dateutil.parser.parse(date+'T'+config_in['night_start'])
        hours_ticks = [time_start + dateutil.relativedelta.relativedelta(hours=i) for i in
                      range(0, int(config_in['night_duration']))]


        plt.subplots_adjust(bottom=0.2)
        plt.xticks(hours_ticks, rotation=80)
        #plt.xticks(minutes_ticks, rotation=80)

        #plt.autoscale(tight=True)

        ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
        #ax.xaxis.set_major_formatter(md.DateFormatter('%d/%m  %H:%M'))
        #ax.xaxis.set_minor_formatter(md.DateFormatter('%M'))

        #ax.xaxis.set_major_locator(md.MinuteLocator(byminute=[0,30], interval=1, tz=None))
        #ax.xaxis.set_minor_locator(md.MinuteLocator(byminute=[10,20,40,50], interval=1, tz=None))
        ax.xaxis.set_major_locator(md.MinuteLocator(byminute=[0,60], interval=1, tz=None))
        ax.xaxis.set_minor_locator(md.MinuteLocator(byminute=[15,30,45], interval=1, tz=None))


        ax.set_title(date)
        ax.set_xlabel('Time (UT)')
        ax.set_ylabel('RV [m/s]')
        ax.plot(plot_bjd, RV, color='black', label='Planet '+planet_name)
        if multiplanets:
            ii = 0
            for pp_name, pp in planets.items():
                if planet_name == pp_name:
                    continue
                ax.plot(plot_bjd, RV_dict[pp_name], color='C'+repr(ii),  label=pp_name)
                ii += 1
            ax.plot(plot_bjd, RV_all, color='red',  label='All signals')

        ax.set_xlim(plot_bjd[0], plot_bjd[-1])
        plt.legend(facecolor='white')

        fig.savefig(yaml_name + '_' + date + '_' + planet_name + '.pdf', bbox_inches='tight', dpi=300)
        plt.close()
