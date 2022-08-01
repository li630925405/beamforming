import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate

fig = plt.figure()

ax = fig.add_subplot(211, projection='polar')
ax2 = fig.add_subplot(212)
ax.set_thetamin(0.0)
ax.set_thetamax(180.0)
# ax.set_rlim(0.0, 10.0)

def tf(source, mic, f):
    c = 343
    r = np.sqrt((source[0] - mic[0]) ** 2 + (source[1] - mic[1]) ** 2)
    k = 2 * np.pi * f / c
    return 1 / (4 * np.pi * r) * np.exp(-1j * k[i] * r[n][m])

# source_position = np.array([[0, 0], [1, 0], [2, 0], ])
N = 10  # number of speakers in array
f = 2000  # frequency to plot
d = 0.08  # distance between each sensor
c = 343  # sound velocity

method = input("method: ")
# method = "LS"

if method == "DS":  
    """DS (delay and sum)"""
    theta = np.deg2rad(30)  # incident angle
    h = []
    # for phi in angle:
    #     A = np.abs(np.sin(N * np.pi * f * d * (np.cos(phi) - np.cos(theta)) / c)
    #                 / (N * np.sin(np.pi * f * d * (np.cos(phi) - np.cos(theta)) / c)))
    #     amplitude.append(10 * np.log10(A))

    for n in range(N):
        h.append(np.exp(1j * 2 * np.pi * n * f * d * np.cos(theta) / c))

    name = "N" + str(N) + "_f" + str(f) + "_d" + str(d) + "_theta" + str(round(np.rad2deg(theta))) + "_" + method + ".jpg"


elif method == 'LS':  
    """LS (least square)"""
    phi1 = np.deg2rad(25)
    phi2 = np.deg2rad(35)
    # phi3 = np.deg2rad(120)
    # phi4 = np.deg2rad(125)
    _d = 2 * np.pi * f * d / c
    Q = np.zeros((N, N))
    p = np.zeros((N))
    
    def calculate_Q(phi):
        s = []
        sH = []
        for i in range(N):
            s.append(np.exp(-1j * i * _d * np.cos(phi)))
            sH.append(np.exp(1j * i * _d * np.cos(phi)))
        s = np.array(s).reshape(N, 1)
        sH = np.array(sH).reshape(1, N)
        return np.matmul(s, sH)
    
    def element_Q(phi, m, n):
        return calculate_Q(phi)[m, n]

    for m in range(N):
        for n in range(N):
            Q[m, n] = integrate.quad(element_Q, 0, np.pi, args=(m, n))[0]

    def calculate_p(phi):
        p = []
        for i in range(N):
            p.append(np.cos(_d * i * np.cos(phi)))
        return p

    def element_p(phi, n):
        return calculate_p(phi)[n]

    for n in range(N):
        p[n] = integrate.quad(element_p, phi1, phi2, args=(n))[0] 
            # + integrate.quad(element_p, phi3, phi4, args=(n))[0]

    h = np.matmul(np.linalg.inv(Q), p)

    name = "N" + str(N) + "_f" + str(f) + "_d" + str(d) + "_phi1" + str(round(np.rad2deg(phi1))) + "_phi2" + str(round(np.rad2deg(phi2)))+ "_" + method + ".jpg"

elif method == "d":
    """discrete version of LS"""
    phi1 = np.deg2rad(25)
    phi2 = np.deg2rad(35)
    # phi3 = np.deg2rad(120)
    # phi4 = np.deg2rad(125)
    _d = 2 * np.pi * f * d / c
    ang = np.deg2rad(np.arange(0, 180, 0.1))
    phi_n = int(180 / 0.1)
    s = np.zeros((N, phi_n))
    p = np.zeros((N, phi_n))

    for a in range(len(ang)):
        for i in range(N):
            s[i][a] = np.exp(-1j * i * _d * np.cos(ang[a]))
            p[i][a] = np.cos(_d * i * np.cos(ang[a]))
    
    s = np.sum(s, axis=1).reshape(N, 1)
    p = np.sum(p, axis=1)
    s_h = s.conjugate().T
    beta = 10e-3
    Q = s @ s_h
    # numpy.linalg.LinAlgError: Singular matrix 
    print(Q.shape)
    print(p.shape)
    h = np.linalg.solve(s @ s_h + beta * np.identity(s.shape[1]), p)
    # h = np.linalg.solve(s @ s_h, p)

else:
    print("no such method")

'''calculate directional response'''
angle = np.deg2rad(np.arange(0, 180, 0.1))
amplitude = []
for phi in angle:
    S = 0
    for n in range(0, N):
        S += h[n] * np.exp(-1j * 2 * np.pi * n * f * d * (np.cos(phi)) / c)
    amplitude.append(S)
amplitude = amplitude / max(amplitude)
amplitude = [10 * np.log10(x) for x in amplitude]


plt.title(name)
ax.plot(angle, amplitude)
ax2.plot(angle, amplitude)

plt.savefig(name, dpi=300)
plt.show()