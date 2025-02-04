from pylab import *

import glob as gl


def readField(time, name, dimens=[-1]):
    p = open(time+'/'+name,'r')
    field = p.read()
    if(dimens==[-1]):
        field = field[field.find('('):field.find('\n)')]
    else:
        field = field[field.find('(')+2:field.find('\n)')]
    field = field.replace('\n',' ')
    field = field.replace('(',' ')
    field = field.replace(';',' ')
    field = field.replace(')',' ')
    field = field.replace('  ',' ')
    field = field.replace('  ',' ')
    if(dimens==[-1]):
        field = array(field.split(' ')[1:])
    else:
        field = array(field.split(' ')[1:-1])
    field = double(field)
    p.close()
    field.shape = dimens
    return(field)


    




procesadores = gl.glob('processor*')


processor_results = []
#[[1,1,''], [1,100,'^'], [100,1,'o'], [100,100,'>']]
studys = [[1,1,'']] 
for study_case in studys:
    processor_results = []
    Ma = study_case[0]
    Tr = 4.0
    Theta = study_case[1]
    rootFolder = 'Ma{}Theta{}_large'.format(Ma, Theta)
    Fv = 0
    Ft = 0
    Fma = 0
    Fdilatacional = []
    Fextensional = []
    print('Caso: {}'.format(study_case))
    times = gl.glob(rootFolder+'/[0-9]*')
    Gamma = []
    vx_dx = []
    vy_dx = []
    vx_dx = []
    vy_dy = []
    nx = []
    ny = []
    Forces_v = []
    Forces_ma = []
    Forces_t = []
    Forces_dilatacional = []
    Forces_extensional = []
    totalS = 0.0
    for time in times:
        try:
            print(time)
            Gamma = array(readField(time, 'magGamma'))
            Sf = readField(time, 'sf',[-1,3])
            magSf = readField(time, 'magsf')
            gradU = readField(time, 'gradU',[-1,9])
            vx_dx = array(gradU[:,0])
            vx_dy = array(gradU[:,1])
            vy_dx = array(gradU[:,3])
            vy_dy = array(gradU[:,4])
            
            Sf = Sf/magSf.reshape([-1,1])
            
            Fma = sum(Ma*Sf[:,1]*Gamma*array(magSf))
            Fv = sum((Sf[:,1]*(2*Tr*vy_dy
                           + (Theta-Tr)*(vx_dx+vy_dy))
                  +Sf[:,0]*(vy_dx+vx_dy)
            )*array(magSf))
            Fdilatacional = sum((Sf[:,1]*(2*Tr*vy_dy))*array(magSf))
            Fextensional = sum((Sf[:,1]*((Theta-Tr)*(vx_dx+vy_dy)))*array(magSf))
            Ft = Fma+Fv
            totalS = sum(magSf)
            Forces_v.append(Fv/totalS)
            Forces_dilatacional.append(Fdilatacional/totalS)
            Forces_extensional.append(Fextensional/totalS)
            Forces_t.append(Ft/totalS)
            Forces_ma.append(Fma/totalS)
            print('Total force: {}     Viscous force:{}'.format(Ft,Fv))
        except Exception as e:
            print(e)
            totalS = 1.0
            Forces_v.append(Fv/totalS)
            Forces_t.append(Ft/totalS)
            Forces_ma.append(Fma/totalS)
            Forces_dilatacional.append(Fdilatacional/totalS)
            Forces_extensional.append(Fextensional/totalS)
            print(e)
    processor_results.append(c_[Forces_t, Forces_v,Forces_ma, Forces_dilatacional, Forces_extensional])
    times = [time.replace(rootFolder+'/','') for time in times]
    processor_results =  sum(processor_results, axis=0)
    times = double(times)
    aa = argsort(times)
    processor_results = processor_results[aa,:]
    times = times[aa]

    plot(times[1:], processor_results[1:,0],'{}k'.format(study_case[2]), label = 'Ftotal')
    plot(times[1:], processor_results[1:,1],'--{}r'.format(study_case[2]), label = 'F_V')
    plot(times[1:], processor_results[1:,2],':{}b'.format(study_case[2]), label = 'F_Ma')
    plot(times[1:], processor_results[1:,3],'-.{}g'.format(study_case[2]), label = 'F_dilatacional')
    plot(times[1:], processor_results[1:,4],'-o{}b'.format(study_case[2]), label = 'F_extensional')

    title('Forces')
    ylabel('F_i')
    xlabel(r'$t\times \omega$')
    legend()


savefig('forces.png')
show()
# p = open(time+'/magGamma','r')
        # field = p.read()
        # field = field[field.find('(')+2:field.find('\n)')]
        # field = field.replace('\n',' ')
        # field = field.replace('(',' ')
        # field = field.replace(';',' ')
        # field = field.replace(')',' ')
        # field = field.replace('  ',' ')
        # field = field.replace('  ',' ')        
        # field = array(field.split(' ')[1:-1])
        # field = double(field)
        # p.close()
